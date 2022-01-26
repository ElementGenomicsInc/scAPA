cd /data

#wget ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
#gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz


#wget ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz
#gunzip Homo_sapiens.GRCh37.87.gtf.gz


#cellranger mkgtf Homo_sapiens.GRCh37.87.gtf Homo_sapiens.GRCh37.87.filtered.gtf \
#                 --attribute=gene_biotype:protein_coding \
#                 --attribute=gene_biotype:lincRNA \
#                 --attribute=gene_biotype:antisense


cellranger mkref --genome=hg19 \
                 --fasta=Homo_sapiens.GRCh37.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh37.87.filtered.gtf \
                 --ref-version=3.0.0 \
                 --nthreads=15