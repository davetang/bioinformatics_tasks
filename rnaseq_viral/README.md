# Mapping bulk RNA-seq reads to the virus genome

Discussion on [mapping to the virus genome](https://github.com/alexdobin/STAR/issues/835) using STAR align:

* If the sample contains both viral and host reads, mapping to the viral genome only may indeed be very slow, as most of the reads will not map.
* It is recommended to map to a combined host/viral genome. Or you can try to reduce the `--seedPerWindowNmax` from the default 50 to 30 or 20 or 10 - please see the discussion here: <https://groups.google.com/d/msg/rna-star/7ZUTnk8_bEI/BKFobC46CgAJ>
* There is no need to concatenate FASTA files; you can list multiple fasta files in the `--genomeFastaFiles` command. You just need to make sure that the chromosome names are distinct in two species.
* Furthermore, it is OK to have only human annotations.
* After getting the sorted BAM and indexing it, use `samtools` to extract the entire viral chromosome:

```console
samtools view -b Aligned.sorted.bam _ViralChr_   >   Aligned.viral.bam
```

* Another approach is to first map to mouse and then using the unmapped mouse reads to map to the viral genome.
* As noted in the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) in `2 Generating genome indexes.`:

> `--sjdbGTFfile` specifies the path to the file with annotated transcripts in the standard GTF format. STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are available. Starting from 2.4.1a, the annotations can also be included on the fly at the mapping step.

* From [STAR run without gtf file](https://github.com/alexdobin/STAR/issues/1455) make sure that `--quantMode TranscriptomeSAM` and/or `GeneCounts` is not used because you need to supply the GTF file that specifies gene/transcript annotations.
* [2pass mode without genome annotation](https://github.com/alexdobin/STAR/issues/1207)

## Data

Get download links for SRR953479.

```console
ffq SRR953479 | grep url
```
```
"urltype": "ftp",
"url": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR953/SRR953479/SRR953479.fastq.gz"
"urltype": "aws",
"url": "s3://sra-pub-src-6/SRR953479/MCMV-infected%20fibroblasts%20110317_s_8_1_fastq.txt.gz"
"urltype": "aws",
"url": "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR953479/SRR953479"
"urltype": "aws",
"url": "s3://sra-pub-zq-8/SRR953479/SRR953479.sralite.1"
"urltype": "gcp",
"url": "gs://sra-pub-zq-106/SRR953479/SRR953479.noqual.1"
"urltype": "ncbi",
"url": "https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-11/SRR000/953/SRR953479/SRR953479.sralite.1"
```

Download RNA-seq data from EBI.

```console
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR953/SRR953479/SRR953479.fastq.gz
```

Download viral reference.

```console
./datasets download virus genome accession NC_004065.1 --include genome cds annotation
unzip -d NC_004065.1 ncbi_dataset.zip
rm ncbi_dataset.zip
```

Download human reference.

```console
datasets download genome accession GCF_000001405.40 --include gff3,rna,cds,protein,genome,seq-report
unzip -d GCF_000001405.40 ncbi_dataset.zip
rm ncbi_dataset.zip
cd GCF_000001405.40

# failed twice
md5sum -c md5sum.txt
```
```
ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna: FAILED
ncbi_dataset/data/GCF_000001405.40/genomic.gff: FAILED
```

## STAR reference

Generate reference without the viral genome.

```console
mkdir -p GRCh38
STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir GRCh38 \
    --genomeFastaFiles /data/genome/GCF_000001405.40/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna \
    --sjdbGTFfile /data/genome/GCF_000001405.40/ncbi_dataset/data/GCF_000001405.40/genomic.gff \
    --sjdbOverhang 100
```

Generate reference with the viral genome.

```console
mkdir -p GRCh38_NC_004065
STAR \
    --runThreadN 4 \
    --runMode genomeGenerate \
    --genomeDir GRCh38_NC_004065 \
    --genomeFastaFiles /data/genome/GCF_000001405.40/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna NC_004065.1/ncbi_dataset/data/genomic.fna \
    --sjdbGTFfile /data/genome/GCF_000001405.40/ncbi_dataset/data/GCF_000001405.40/genomic.gff \
    --sjdbOverhang 100
```

## STAR align

Map to single reference.

```console
STAR \
   --genomeDir /data/index/GRCh38 \
   --runThreadN 4 \
   --readFilesCommand gunzip -c \
   --readFilesIn SRR953479.fastq.gz \
   --outSAMtype BAM SortedByCoordinate \
   --twopassMode Basic \
   --outFileNamePrefix SRR953479.star.single.
```

Map to dual reference.

```console
STAR \
   --genomeDir /data/index/GRCh38_NC_004065 \
   --runThreadN 4 \
   --readFilesCommand gunzip -c \
   --readFilesIn SRR953479.fastq.gz \
   --outSAMtype BAM SortedByCoordinate \
   --twopassMode Basic \
   --outFileNamePrefix SRR953479.star.dual.
```

