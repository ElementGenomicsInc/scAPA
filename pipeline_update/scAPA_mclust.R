require("scAPA")

script.args = commandArgs(trailingOnly = TRUE)
path.to.files = script.args[1]
c = as.integer(script.args[2])
org = script.args[3]

setwd(path.to.files)

peaks.wig <- read.delim(file = "intersected.wig", header = F)
peaks.wig <- split(x = peaks.wig, f = peaks.wig$V10, drop = T)

if(c > 1){
  bed <- plyr::rbind.fill(parallel::mclapply(1:length(peaks.wig),
                                             FUN = creat_mclus,
                                             mc.cores = c,
                                             mc.preschedule = T))
}
if(c == 1){
  bed <- plyr::rbind.fill(lapply(X = 1:length(peaks.wig),
                                 FUN = creat_mclus))
}