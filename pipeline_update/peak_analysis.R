require("scAPA")

script.args = commandArgs(trailingOnly = TRUE)
path.to.files = script.args[1]
bedtools.path = script.args[2]
org = script.args[3]

setwd(paste(path.to.files, "temp", sep="/"))
merge_peaks(bedtools.path = bedtools.path, peaks.file = "Peakfile", path = "./")
peaks.bed <- intersect_peaks(org = org, bed.name = "./merge.peakfile.bed",
                             path = "", bedtools.path = bedtools.path,
                             introns = F)
write.bed(.x = peaks.bed, f = "./peaks.bed")