rsem-prepare-reference \
-p 3 \
--gtf  $gtf \
$genome \
$out

STAR --genomeDir $genome \
--readFilesIn $r1 $r2  \
--outFileNamePrefix  ${outdir}/$name \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 5 --limitOutSJcollapsed 5000000 \
--quantMode TranscriptomeSAM GeneCounts

rsem-calculate-expression --paired-end --no-bam-output --alignments -p 5 \
-q $bam \
$genome  ${outdir}/${name}