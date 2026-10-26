"""Offline regression tests. Run: python3 -m unittest discover -s tests -v."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "figshare-dl"

# Replace only the network boundary; the real Bash script, jq and hash tools run.
FAKE_CURL = r'''#!/usr/bin/env python3
import json
import os
from pathlib import Path
import sys
from urllib.parse import parse_qs, urlsplit

args = sys.argv[1:]
fixture = json.loads(Path(os.environ["FIXTURE"]).read_text())
with open(os.environ["CURL_LOG"], "a") as log:
    log.write(json.dumps(args) + "\n")
url = next(arg for arg in args if arg.startswith("https://"))
path = urlsplit(url).path
if "-o" in args:
    entry = next(f for f in fixture["files"] if f["download_url"] == url)
    output = Path(args[args.index("-o") + 1])
    data = entry.get("body", "hello").encode()
    offset = output.stat().st_size if output.exists() else 0
    with output.open("ab") as target:
        if entry.get("fail_download"):
            target.write(data[offset:offset + 2])
            sys.exit(22)
        target.write(data[offset:])
    sys.exit(0)
if path.endswith("/search"):
    json.loads(args[args.index("-d") + 1])
    result = [{"id": 123}]
elif not path.startswith("/v2/articles/123/"):
    sys.exit(22)
elif path.endswith("/versions"):
    result = fixture.get("versions", [{"version": 1}, {"version": 2}])
elif path.endswith("/files"):
    page = int(parse_qs(urlsplit(url).query)["page"][0])
    if page == fixture.get("fail_page"):
        sys.exit(22)
    if fixture.get("malformed_files"):
        print("not JSON")
        sys.exit(0)
    if "pages" in fixture:
        result = fixture["pages"][page - 1]
    else:
        result = fixture["files"]
else:
    result = {"title": "Test article", "doi": "10.6084/m9.figshare.123"}
print(json.dumps(result))
'''


class FigshareDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        self.output = self.root / "output"
        self.output.mkdir()
        self.fixture_path = self.root / "fixture.json"
        self.log = self.root / "curl.jsonl"
        self.fixture = {"files": [self.file()]}
        curl = self.bin / "curl"
        curl.write_text(FAKE_CURL)
        curl.chmod(0o755)
        self.env = dict(os.environ, PATH=f"{self.bin}:{os.environ['PATH']}",
                        FIXTURE=str(self.fixture_path), CURL_LOG=str(self.log))

    @staticmethod
    def file(name="data.txt", body="hello", **overrides):
        entry = {"name": name, "body": body, "size": len(body.encode()),
                 "computed_md5": hashlib.md5(body.encode()).hexdigest(),
                 "download_url": "https://download.test/" + str(len(name))}
        entry.update(overrides)
        return entry

    def run_script(self, article="123", options=(), expected=0):
        self.fixture_path.write_text(json.dumps(self.fixture))
        result = subprocess.run(["bash", str(SCRIPT), *options, article, str(self.output)],
                                env=self.env, cwd=self.root, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20)
        self.assertEqual(result.returncode, expected, result.stdout + result.stderr)
        return result

    def downloads(self):
        return [json.loads(line) for line in self.log.read_text().splitlines()
                if "-o" in json.loads(line)]

    def test_download_then_skip_verified_file(self):
        self.run_script()
        self.assertEqual((self.output / "data.txt").read_text(), "hello")
        self.assertIn("have", self.run_script().stdout)
        self.assertEqual(len(self.downloads()), 1)

    def test_resume_after_download_failure(self):
        self.fixture["files"][0]["fail_download"] = True
        self.run_script(expected=1)
        self.assertEqual((self.output / "data.txt.part").read_text(), "he")
        self.fixture["files"][0]["fail_download"] = False
        self.run_script()
        self.assertEqual((self.output / "data.txt").read_text(), "hello")
        self.assertFalse((self.output / "data.txt.part").exists())

    def test_complete_partial_is_promoted_without_downloading(self):
        (self.output / "data.txt.part").write_text("hello")
        self.run_script()
        self.assertEqual(self.downloads(), [])
        self.assertEqual((self.output / "data.txt").read_text(), "hello")

    def test_missing_checksums_still_require_correct_size(self):
        for checksum in (None, ""):
            with self.subTest(checksum=checksum):
                self.fixture["files"][0]["computed_md5"] = checksum
                (self.output / "data.txt").write_text("h")
                self.run_script()
                self.assertEqual((self.output / "data.txt").read_text(), "hello")
        self.assertEqual(len(self.downloads()), 2)

    def test_truncated_download_without_checksum_is_rejected(self):
        self.fixture["files"][0].update(computed_md5=None, body="h")
        self.run_script(expected=1)
        self.assertFalse((self.output / "data.txt").exists())

    def test_failed_replacement_preserves_old_file(self):
        (self.output / "data.txt").write_text("old contents")
        self.fixture["files"][0]["body"] = "wrong"
        self.run_script(expected=1)
        self.assertEqual((self.output / "data.txt").read_text(), "old contents")
        self.assertFalse((self.output / "data.txt.part").exists())

    def test_bad_complete_partial_is_restarted(self):
        (self.output / "data.txt.part").write_text("wrong")
        self.run_script()
        self.assertEqual((self.output / "data.txt").read_text(), "hello")

    def test_unusual_filenames_and_uppercase_checksum(self):
        name = "-file\twith\\slashes\nand newline\n"
        self.fixture["files"] = [self.file(name)]
        self.fixture["files"][0]["computed_md5"] = self.fixture["files"][0]["computed_md5"].upper()
        self.run_script()
        self.assertEqual((self.output / name).read_text(), "hello")
        self.run_script()
        self.assertEqual(len(self.downloads()), 1)

    def test_file_list_failure_is_reported(self):
        self.fixture["fail_page"] = 1
        self.run_script(expected=1)
        self.assertEqual(self.downloads(), [])

    def test_failure_on_later_page_prevents_any_downloads(self):
        self.fixture.update(pages=[[self.file()] * 1000], fail_page=2)
        self.run_script(expected=1)
        self.assertEqual(self.downloads(), [])

    def test_invalid_file_list_is_reported(self):
        self.fixture["malformed_files"] = True
        self.run_script(expected=1)
        self.assertEqual(self.downloads(), [])

    def test_empty_file_list_succeeds(self):
        self.fixture["files"] = []
        self.run_script()
        self.assertEqual(self.downloads(), [])

    def test_empty_version_list_has_clear_error(self):
        self.fixture["versions"] = []
        self.assertIn("no public versions", self.run_script(expected=1).stderr)

    def test_urls_and_version_selection(self):
        for article in ("https://figshare.com/articles/dataset/2024/123/1?file=9",
                        "https://figshare.com/articles/123/1",
                        "https://api.figshare.com/v2/articles/123/versions/1",
                        "https://figshare.com/articles/dataset/slug/123/1/",
                        "https://doi.org/10.6084/m9.figshare.123.v1"):
            with self.subTest(article=article):
                self.assertIn("version 1", self.run_script(article, options=("-v", "01")).stdout)
        self.assertIn("version 2", self.run_script("https://figshare.com/articles/dataset/2024/123").stdout)

    def test_version_conflict_and_invalid_numbers(self):
        self.run_script("https://figshare.com/articles/dataset/slug/123/1", options=("-v2",), expected=2)
        for version in ("", "abc", "0"):
            with self.subTest(version=version):
                self.run_script(options=("-v", version), expected=2)
        self.run_script(options=("-v", "999999999999999999999999999"), expected=1)

    def test_doi_is_json_escaped(self):
        self.run_script('10.1234/a"b')
        calls = [json.loads(line) for line in self.log.read_text().splitlines()]
        search = next(call for call in calls if "-d" in call)
        self.assertEqual(json.loads(search[search.index("-d") + 1]), {"doi": '10.1234/a"b'})

    def test_unsafe_and_conflicting_filenames_are_rejected(self):
        for names in (("..",), (".",), ("path/",), ("bad\x00name",),
                      ("x", "folder/x"), ("x", "x.part")):
            with self.subTest(names=names):
                self.fixture["files"] = [self.file(name) for name in names]
                self.run_script(expected=1)
        self.assertEqual(self.downloads(), [])

    def test_existing_directory_is_not_treated_as_file(self):
        (self.output / "data.txt").mkdir()
        self.run_script(expected=1)
        self.assertEqual(list((self.output / "data.txt").iterdir()), [])

    def test_md5_fallback(self):
        # Isolate PATH to exercise the macOS command branch even on Linux.
        for command in ("bash", "python3", "jq", "mkdir", "mktemp", "rm", "wc", "mv"):
            (self.bin / command).symlink_to(shutil.which(command))
        md5 = self.bin / "md5"
        md5.write_text("#!/usr/bin/env python3\nimport hashlib, sys\n"
                       "assert sys.argv[1:] == ['-q']\n"
                       "print(hashlib.md5(sys.stdin.buffer.read()).hexdigest())\n")
        md5.chmod(0o755)
        self.env["PATH"] = str(self.bin)
        self.run_script()
        self.assertEqual((self.output / "data.txt").read_text(), "hello")


if __name__ == "__main__":
    unittest.main()
