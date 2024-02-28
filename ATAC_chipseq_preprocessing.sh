
# Align to genome
align_pe() { genome=/mnt/3/ywlai_genome/genome/mm10/mm10_bowtie2_index/mm10; r1=$1; r2=`echo $r1 | sed 's/_1./_2./g'`
   out=`echo $r1 | sed 's/.fastq//g'`; bowtie2 --very-sensitive -p 2 --no-unal -x $genome -1 $r1 -2 $r2 -S $out".bowtie2.sam" 
}
  export -f align_pe


# filter, sorting and generating CPM normalized bigwig files
samtools view -q 30 -bS -F 0x04 $f | samtools sort -@ 2 > $g
bamCoverage -p 20 --bam f -o file.bw --binSize 10 --normalizeUsing CPM

# peak calling for ATAC-seq
macs2 callpeak -B --nomodel --keep-dup 1 -g mm --call-summits -t t -f BAM --outdir out -n name -q 0.01

# peak calling for ChIP-seq
macs2 callpeak --nomodel  -B --keep-dup 1 -g mm --call-summits -t t -c c -f BAM \
--outdir outdir -n name -q 0.01

# Differential peaks comparing H3K27ac ChIP-seq at 40 hours post-PHx and time 0 
macs2 bdgdiff --t1 $t1 --t2 $t2 --c1 $c1 --c2 $c2 --outdir outdir --o-prefix prefix
