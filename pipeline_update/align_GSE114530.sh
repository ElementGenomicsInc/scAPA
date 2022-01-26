mkdir -p /data/GSE114530
cd /data/GSE114530

#fastq-dump --split-files SRR7171583

#for read in 1 2
#do
#  mv SRR7171583_$read.fastq SRR7171583_S1_L001_R${read}_001.fastq
#done

cellranger count --id SRR7171583 --transcriptome /data/hg19 --fastqs /data/GSE114530 --expect-cells 7000 --localcores 15