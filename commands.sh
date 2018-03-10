python /home/ubuntu//Pipeline/demultiplex.py --fastq1 raw_reads/160811Gil_D16-8293_2_sequence.fastq.gz --fastq2 raw_reads/160811Gil_D16-8293_1_sequence.fastq.gz --prefix Trm7_Li_B/intermediates/Trm7_Li_B --metrics Trm7_Li_B/logs/demux.metrics.txt

fastqc Trm7_Li_B/intermediates/Trm7_Li_B.r1.fastq.gz -q -o Trm7_Li_B/fastqc/ >> Trm7_Li_B/logs/Trm7_Li_B.r1.fastqc.log

fastqc Trm7_Li_B/intermediates/Trm7_Li_B.r2.fastq.gz -q -o Trm7_Li_B/fastqc/ >> Trm7_Li_B/logs/Trm7_Li_B.r2.fastqc.log

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 --quality-cutoff 6 -m 18 -a NNNNNGATCGTCGGACTGTAGAACTCTGAACGTGTAG -g GCACCCGAGAATTCCA -g GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA -g CAAGCAGAAGACGGCATACGAGATAACTTG -a TGGAATTCTCGGG -a CCAAGGAACTCCAGTCAC -a TGGAATTCTCGGGTGCCAAGGCACTCCAGTCACCAAGTTATCTCGTATGCCGTCTTCTGCTTG -o Trm7_Li_B/intermediates/Trm7_Li_B.r1.trimmed.fastq.gz -p Trm7_Li_B/intermediates/Trm7_Li_B.r2.trimmed.fastq.gz Trm7_Li_B/intermediates/Trm7_Li_B.r1.fastq.gz Trm7_Li_B/intermediates/Trm7_Li_B.r2.fastq.gz >> Trm7_Li_B/logs/cutadapt.log

cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 5 --quality-cutoff 6 -m 18 -a TGGAATTCTCGGG -a CCAAGGAACTCCAGTCAC -a TGGAATTCTCGGGTGCCAAGGCACTCCAGTCACCAAGTTATCTCGTATGCCGTCTTCTGCTTG -o Trm7_Li_B/intermediates/Trm7_Li_B.r1.trimmed.round2.fastq.gz -p Trm7_Li_B/intermediates/Trm7_Li_B.r2.trimmed.round2.fastq.gz Trm7_Li_B/intermediates/Trm7_Li_B.r1.trimmed.fastq.gz Trm7_Li_B/intermediates/Trm7_Li_B.r2.trimmed.fastq.gz >> Trm7_Li_B/logs/cutadapt.log

fastqc Trm7_Li_B/intermediates/Trm7_Li_B.r1.trimmed.round2.fastq.gz -q -o Trm7_Li_B/fastqc/ >> Trm7_Li_B/logs/Trm7_Li_B.r1.trimmed.round2.fastqc.log

fastqc Trm7_Li_B/intermediates/Trm7_Li_B.r2.trimmed.round2.fastq.gz -q -o Trm7_Li_B/fastqc/ >> Trm7_Li_B/logs/Trm7_Li_B.r2.trimmed.round2.fastqc.log

STAR --runMode alignReads --runThreadN 2 --outFilterMismatchNoverLmax 0.09 --genomeDir /home/ubuntu//Pipeline/sacCer3/sacCer3_STAR_index/ --genomeLoad NoSharedMemory --readFilesIn Trm7_Li_B/intermediates/Trm7_Li_B.r1.trimmed.round2.fastq.gz Trm7_Li_B/intermediates/Trm7_Li_B.r2.trimmed.round2.fastq.gz --readFilesCommand zcat --outSAMunmapped Within --outFilterMultimapNmax 20 --outFilterMultimapScoreRange 1 --outFileNamePrefix Trm7_Li_B/Trm7_Li_B.genome_mapped.bam --outSAMattributes All --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outFilterType BySJout --outReadsUnmapped Fastx --outFilterScoreMin 10 --outSAMattrRGline ID:foo --alignEndsType EndToEnd > Trm7_Li_B/Trm7_Li_B.genome_mapped.bam

samtools view -q 200 -hb Trm7_Li_B/Trm7_Li_B.genome_mapped.bam > Trm7_Li_B/Trm7_Li_B.genome_mapped.unique.bam

python /home/ubuntu/Pipeline/gscripts-1.1/gscripts/clipseq/barcode_collapse_pe.py --bam Trm7_Li_B/Trm7_Li_B.genome_mapped.unique.bam --out_file Trm7_Li_B/mapped/genome_mapped.rmDup.bam --metrics_file Trm7_Li_B/logs/rmDup.log

java -jar /home/ubuntu//picard/picard.jar SortSam INPUT=Trm7_Li_B/mapped/genome_mapped.rmDup.bam OUTPUT=Trm7_Li_B/Trm7_Li_B.rmDup.bam SO=coordinate VALIDATION_STRINGENCY=SILENT > Trm7_Li_B/logs/SortSam.log

samtools index Trm7_Li_B/Trm7_Li_B.rmDup.bam

samtools view -hb -f 128 Trm7_Li_B/Trm7_Li_B.rmDup.bam > Trm7_Li_B/Trm7_Li_B.rmDup.r2.bam

clipper -b Trm7_Li_B/Trm7_Li_B.rmDup.r2.bam -s sacCer3-o Trm7_Li_B/Trm7_Li_B.rmDup.r2.peaks.bed

python /home/ubuntu/Pipeline/gscripts-1.1/gscripts//clipseq/fix_scores.py --bed Trm7_Li_B/Trm7_Li_B.rmDup.r2.peaks.bed --out_file Trm7_Li_B/Trm7_Li_B.peaks.fixed.bed

