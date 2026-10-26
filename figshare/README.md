# figshare-dl

`figshare-dl` is a Bash script that downloads the files belonging to one public Figshare article version. It accepts an article ID, an article URL, or a DOI, then checks each downloaded file against the size and MD5 checksum reported by Figshare.

Downloads can be resumed by running the same command again. Files that already pass verification are skipped. The script uses public API endpoints and does not require a Figshare account or API token.

## Requirements

- **Bash**: run the script directly or with `bash`; do not invoke it with `sh`.
- **curl**: makes API requests and downloads files.
- **jq**: reads and validates JSON responses.
- **md5sum** or **md5**: computes file checksums. The script prefers `md5sum` and falls back to `md5`.
- Standard command-line utilities: `mkdir`, `mktemp`, `rm`, `mv`, `wc`, and `cut`.
- Internet access and enough disk space for the downloads. Replacing an existing file temporarily requires space for both the old file and its replacement.

Python 3 is needed only to run the regression tests.

## Quick start

From the directory containing the script:

```bash
chmod +x figshare-dl
./figshare-dl -v 3 5616409
```

This downloads version 3 of article `5616409` into `figshare_5616409_v3/`, relative to the current directory.

To use an explicit output directory:

```bash
./figshare-dl -v 3 5616409 ./downloads
```

## Usage

```text
figshare-dl [-v N] <article-id | figshare URL | DOI> [outdir]
```

| Argument or option | Meaning |
| --- | --- |
| Article reference | One numeric article ID, article URL, bare DOI, or DOI URL. |
| `outdir` | Optional destination directory; created if necessary. Defaults to `figshare_<id>_v<version>`. |
| `-v N` or `-vN` | Select a positive integer version number. |
| `-h`, `--help` | Print the usage line and exit. |
| `--` | End option parsing; subsequent arguments are positional. |

Examples:

```bash
# Download the latest public version using an article ID.
./figshare-dl 5616409

# Resolve an article using its DOI.
./figshare-dl '10.6084/m9.figshare.5616409'

# Select a version using the DOI suffix.
./figshare-dl 'https://doi.org/10.6084/m9.figshare.5616409.v3'

# Use a version-specific API URL as the article reference.
./figshare-dl 'https://api.figshare.com/v2/articles/5616409/versions/3'

# Use a destination containing spaces.
./figshare-dl -v3 5616409 './downloaded files'
```

Website URLs of the form `https://figshare.com/articles/<type>/<slug>/<id>` are also accepted, with an optional `/<version>` suffix. The script ignores query strings and fragments when parsing article URLs. Quote URLs so shell characters such as `?` and `&` are passed literally.

### Version selection

The script uses `-v` if supplied, otherwise a version embedded in the URL or DOI, otherwise the highest public version returned by the API. An explicit `-v` must agree with any embedded version; a conflict is an error rather than an override. Leading zeros are removed, so `-v 03` selects version 3.

The chosen version must appear in the article's public version list. To repeat a download of the same version later, specify that version: an unversioned reference may select a newer version on a future run.

## Figshare API background

Figshare stores research outputs as items, called **articles** in its API. An article contains descriptive metadata and may have multiple files. The numeric article ID identifies the item; each file has its own ID and download URL. Public metadata and publicly available files can be retrieved without authentication; private account operations require authentication. See Figshare's [API introduction and examples](https://info.figshare.com/user-guide/how-to-use-the-figshare-api/).

Figshare items can have multiple published versions. A standard Figshare base DOI points to the latest public version, while a suffix such as `.v3` identifies a particular version. Institutional or publisher-supplied DOIs may follow different conventions. See [How versioning works](https://info.figshare.com/user-guide/how-versioning-works/).

The script uses the API base URL `https://api.figshare.com/v2`. Metadata requests return JSON, which `jq` extracts and validates. File contents are fetched separately using the returned `download_url`.

### Endpoints used by the script

The following paths are relative to the API base URL. Their definitions are available in the [Figshare API reference](https://docs.figshare.com/v2/) and its [OpenAPI schema](https://docs.figshare.com/docs/versions/swagger_v2.0.json).

| Method | Endpoint | Use in `figshare-dl` |
| --- | --- | --- |
| `POST` | `/articles/search` | Resolve a DOI using a JSON body such as `{"doi":"10.6084/m9.figshare.5616409"}`. The script takes the first result's article ID. |
| `GET` | `/articles/{id}/versions` | Find the available public versions and validate the selected version. |
| `GET` | `/articles/{id}/versions/{version}` | Retrieve the selected version's title and DOI for display. |
| `GET` | `/articles/{id}/versions/{version}/files?page=N&page_size=1000` | Retrieve file metadata in pages. |

The script follows each file's `download_url`, including HTTP redirects, rather than constructing a file URL from its name. The relevant metadata fields are:

| Field | Use |
| --- | --- |
| `download_url` | Source of the file contents. |
| `name` | Local filename, after removing any directory prefix. |
| `size` | Expected byte count. |
| `computed_md5` | Expected MD5 digest, when present and nonempty. |

The script does not fall back to `supplied_md5`. When `computed_md5` is missing or empty, it checks size only; this cannot detect different contents of the same size.

To inspect a file list yourself without downloading its contents:

```bash
curl -fsS \
  'https://api.figshare.com/v2/articles/5616409/versions/3/files?page=1&page_size=1000' \
  | jq .
```

## How the script works

1. **Parse the command line.** Validate the arguments and check that the required tools are available.
2. **Resolve the article and version.** Extract an ID from a numeric reference or article URL, or search for a DOI. Fetch the public versions and select one.
3. **Retrieve all file metadata.** Request pages of up to 1,000 files until a page contains fewer than 1,000 entries. Finish this step before downloading anything, so a failed later page cannot be mistaken for a complete listing.
4. **Validate the file list.** Check URLs, names, sizes and checksum formats. Reject conflicting local filenames, including a real filename that would overlap another file's `.part` path.
5. **Download and verify each file.** Skip verified destinations, resume partial downloads when possible, and move completed files into place only after verification.
6. **Report the result and clean up.** Remove temporary metadata files when the script exits. Retain partial downloads after transfer failures so another run can resume them.

Temporary metadata is stored in a directory created by `mktemp`. The download loop reads fields separated by NUL bytes, which preserves tabs, newlines and backslashes in filenames. Files are downloaded sequentially.

### Existing files and interrupted downloads

For a file named `data.csv`, the final path is `<outdir>/data.csv` and the temporary download path is `<outdir>/data.csv.part`.

| Local state | Action |
| --- | --- |
| Final file matches size and any published checksum | Skip it. |
| Final file fails verification | Keep it until a verified replacement is ready. |
| `.part` file already passes verification | Move it to the final path without downloading again. |
| `.part` file is shorter than expected | Ask curl to resume from its current length. |
| `.part` file is at least the expected size but fails verification | Delete the partial file and restart the download. |
| Transfer fails | Keep the partial file, report failure, and continue with the next file. |
| Completed transfer fails size or MD5 verification | Delete the partial file, report failure, and continue with the next file. |
| Destination or partial path is a directory or symlink | Report failure for that file. |

Downloads use curl's `--retry 3` for errors curl considers retryable. Resuming depends on the download server supporting range requests. If repeated attempts cannot resume a particular file, remove only its `.part` file and rerun to start that file again.

### Output and exit status

The script prints the article title, selected version, DOI and destination, followed by per-file messages:

| Message | Meaning |
| --- | --- |
| `have` | An existing final file passed verification. |
| `get` | A download is starting or resuming. |
| `ok` | A completed partial file passed verification and was moved into place. |
| `FAIL` | That file could not be downloaded or verified. |

A newly completed file without a published checksum also produces a size-only verification notice on stderr.

Exit status `0` means the requested work completed successfully, including an article with no listed files. Argument errors normally return `2`; handled lookup, metadata or download failures return `1`. Commands that fail outside explicit error handling can propagate their own nonzero status. Automation should treat **any nonzero status as failure**.

## Function guide

Each function also has a description directly above its declaration in [figshare-dl](figshare-dl).

| Function | Responsibility |
| --- | --- |
| `fail` | Print an error and exit with the supplied status. |
| `normalize_number` | Validate a positive integer and strip leading zeros without integer arithmetic. |
| `parse_arguments` | Read options, the article reference and the output directory. |
| `check_dependencies` | Check curl and jq, and select the checksum utility. |
| `resolve_article` | Determine the article ID and any embedded version; detect version conflicts. |
| `select_version` | Validate the chosen public version and build its API URL. |
| `list_files` | Fetch every metadata page and emit one JSON object per file. |
| `prepare_file_list` | Validate metadata and build the NUL-separated download list. |
| `file_md5` | Compute a file's MD5 digest using the selected utility. |
| `file_matches` | Check a local file's size and optional checksum. |
| `download_files` | Skip, resume, verify and finalize individual files. |
| `main` | Coordinate the workflow, arrange cleanup and exit with the result. |

## Scope and limitations

- One public article version is processed per invocation. Collections, projects, private items and authenticated downloads are not supported.
- DOI resolution depends on the public search endpoint returning a matching article. The script does not resolve arbitrary DOI redirects itself.
- Files are stored in one output directory. Directory prefixes in API filenames are discarded; names that then collide cause the run to fail before downloading.
- The output directory is not a mirror: unrelated files, or files left over from another version, are not removed. Separate directories per version avoid mixing their contents.
- Externally hosted files must provide a usable download URL and metadata that pass the same checks. The script has no special handling for linked files with unknown sizes.
- API metadata is temporary; the script downloads file contents but does not save a permanent metadata record.

## Tests

Run the offline regression suite from this directory:

```bash
python3 -m unittest discover -s tests -v
```

The tests in [tests/test_figshare_dl.py](tests/test_figshare_dl.py) replace curl with a local stand-in while running the real script, jq and checksum commands. They cover version selection, DOI encoding, pagination failures, unusual filenames, checksum and size checks, interrupted transfers, existing files, and the `md5` fallback branch. They do not require network access and do not verify the current availability of Figshare's service.

For a Bash syntax check:

```bash
bash -n figshare-dl
```
