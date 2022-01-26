path_to_files=$1
f=$2

cd $path_to_files

mkdir -p temp Log.files

sample=${f%.bam}

#~/Drop-seq_tools-2.2.0/FilterBam -TAG_RETAIN UB -I $f -O temp/UB.$sample.bam > Log.files/FilterBAM.$sample.out

#samtools index temp/UB.$sample.bam > Log.files/Index.$sample.out

umi_tools dedup -I temp/UB.$sample.bam -S temp/dedup.$sample.bam --method=unique --extract-umi-method=tag --umi-tag=UB --cell-tag=CB > Log.files/umi_tools.$sample.out

test -s temp/dedup.$sample.bam && (rm temp/UB.$sample.bam)
