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
