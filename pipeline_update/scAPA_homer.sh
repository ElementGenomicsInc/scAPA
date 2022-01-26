path_to_files=$1

cd $path_to_files

~/homer/bin/makeTagDirectory temp/Tagdirectory temp/dedup.*.bam > Log.files/makeTagDirectory.out
test -d temp/Tagdirectory && (rm temp/dedup.*.bam)

~/homer/bin/findPeaks temp/Tagdirectory -size 50 -fragLength 100 -minDist 1 -strand separate -o temp/Peakfile > Log.files/findPeaks.out
