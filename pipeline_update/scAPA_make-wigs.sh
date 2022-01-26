path_to_files=$1

cd $path_to_files

samtools merge merged.bam dedup*.bam > Log.files/samtoolsmerge.out

bedtools genomecov -ibam merged.bam -bg -strand + | awk 'BEGIN {OFS = \"\t\"}{print $1 , $2, $3, $4, \".\", \"+\"}' > merged.wig
bedtools genomecov -ibam merged.bam -bg -strand - | awk 'BEGIN {OFS = \"\t\"}{print $1, $2, $3, $4, \".\", \"-\"}' >> merged.wig
test -s merged.wig && (rm merged.bam)

bedtools intersect -s -wb -b peaks.bed -a merged.wig > intersected.wig
test -s intersected.wig && (rm merged.wig)