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
