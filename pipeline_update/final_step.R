require("scAPA")

script.args = commandArgs(trailingOnly = TRUE)
path.to.files = script.args[1]

setwd(path.to.files)

bam.cluster.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                full.names = F)
cellnames <- gsub(x = bam.cluster.files, pattern = ".bam", replacement = "")
bam.cluster.files <- list.files(path = "./CellBams", pattern = ".bam$",
                                full.names = T)
counts <- Rsubread::featureCounts(files = bam.cluster.files, isGTFAnnotationFile = F,
                                  strandSpecific = 1, annot.ext = utr.saf,
                                  largestOverlap = T, nthreads = c)
co <- cbind.data.frame(rownames(counts$counts), counts$counts)
colnames(co) <- c("Peak_ID", cellnames)
meta <- counts$annotation
meta <- meta[, c(2, 3, 4, 1, 6, 5)]
metadata <- read_down.seq(saf = utr.saf, char.length.path = char.length.path,
                          fasta.path = fasta.path, chr.modify = T)
aseq <- metadata[, c(4, 6)]
a <- set_scAPAList(.cells.counts = co, .row.Data = meta, .down.seq = aseq,
                   .cluster.anot = full.list.cells)
rm(list = c("co", "meta", "aseq"))
saveRDS(object = a, file = "../outs/Peaks.RDS")
  counts_int <- Rsubread::featureCounts(files = bam.cluster.files,
                                        isGTFAnnotationFile = F,
                                        strandSpecific = 1,
                                        annot.ext = int.saf,
                                        largestOverlap = T,
                                        nthreads = c)
  co_int <- cbind.data.frame(rownames(counts_int$counts),
                             counts_int$counts)
  colnames(co_int) <- c("Peak_ID", cellnames)
  meta_int <- counts_int$annotation
  meta_int <- meta_int[, c(2, 3, 4, 1, 6, 5)]
  metadata_int <- read_down.seq(saf = int.saf,
                                char.length.path = char.length.path,
                                fasta.path = fasta.path,
                                chr.modify = T)
  aseq_int <- metadata_int[, c(4, 6)]
  a.int <- set_scAPAList(.cells.counts = co_int, .row.Data = meta_int,
                         .down.seq = aseq_int,
                         .cluster.anot = full.list.cells)
  rm(list = c("co_int", "meta_int", "aseq_int"))
  saveRDS(object = a.int, file = "../outs/Peaks.int.RDS")
  
rm(list = c("counts"))
system(command = "rm -r ./CellBams")


# Peak filtering ----------------------------------------------------------
a <- calc_clusters_counts(a)
a <- calc_cpm(a)
keep.cpm <- rowSums(a@norm$CPM) > CPM.cuttoff
a.fil <- a[keep.cpm, ]

# Internal priming filtering
IP.seq <- paste0(rep("A", times = A.number), collapse = "")
a.fil <- filter_IP(x = a.fil, int.priming.seq = IP.seq,
                   left = filter.border.left, right = filter.border.right)

# Statistical analysis ----------------------------------------------------
results <- set_scAPAresults(a.fil)
# Chi-squre test for APA
results <- test_APA(results, clus = "all")
# For APA with more than one peak and whose p-value is < sig.level,
# perform chi-squered test for goodness of fit
results <- test_peaks(results, clus = "all", sig.level = 0.05)

# Inferring global trends -------------------------------------------------
results <- calc_pAi_mat(results)
results <- calc_p_pui_mat(results)
saveRDS(object = results, file = "../outs/results.RDS")

# Write final outputs -----------------------------------------------------
setwd("../")
results.output <- disply_results(x = results, org = org)
a.fil <- annotate(a.fil, org = org)
  cells_pui <- as.data.frame(results@ppui.cells)
  cells_pui$cell <- gsub(pattern = "_proximal_PUI", replacement = "",
                         x = rownames(cells_pui))
  cells_pui$sample <- gsub(pattern = "_.*", replacement = "",
                           x = cells_pui$cell)
  cells_pui$cell <- gsub(pattern = ".*_", replacement = "",
                         x = cells_pui$cell)
  cells_pui <- cells_pui[, c(3, 2, 1)]
  colnames(cells_pui) <- c("Sample", "Cell_BC", "Mean_Proximal_PUI")
  write.table(x = cells_pui, file = "./outs/Mean.Cell.PPUI.txt",
              quote = F, sep = "\t", col.names = T, row.names = F)

write.table(x = results.output, file = "./outs/APA.events.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.table(x = a.fil@row.Data, file = "./outs/ThreeUTR.peaks.txt",
            quote = F, sep = "\t", col.names = T, row.names = F)

write.peaks(.x = results, .f = "./outs/UTRs.with.multiple.peaks.txt")

file.create("./outs/summary.UTR.txt")
cat("Number of peaks passed fillter:\t", nrow(a.fil@clus.counts),
    "\nNumber of significant (FDR < 5%) APA events:",
    sum(results@pvalues[[1]][, 2] < 0.05, na.rm = T), "\n",
    file = "./outs/summary.UTR.txt")

# Plot --------------------------------------------------------------------
sig.utrs <- as.vector(results@pvalues$all[,2] < 0.05)
sig <- results[which(sig.utrs),]
tidy.ppui <- as.data.frame(sig@ppui.clus)
colnames(tidy.ppui) <- gsub(x = colnames(tidy.ppui),
                            pattern = "_proximal_PUI", replacement = "")
tidy.ppui <- tidyr::gather(data = tidy.ppui)
colnames(tidy.ppui) <- c("Cluster", "value")
pdf("./outs/Proximal.PUI.ECDF.pdf")
p <- ggplot2::ggplot(data = tidy.ppui,
                     ggplot2::aes(x = value, color = Cluster))
p <- p + ggplot2::stat_ecdf(size = 1)
p <- p + ggplot2::theme_bw()
p <- p + ggplot2::xlab("Proximal peak usage index")
p <- p + ggplot2::ylab("Cumulative fraction")
print(p)
dev.off()

# Introns -----------------------------------------------------------------
# Peak filtering ----------------------------------------------------------
 a.int <- calc_clusters_counts(a.int)
  results.int <- set_scAPAresults(x = a.int, int = T, cpm = ICPM.cuttoff,
                                  counts = IC.cuttoff)
  
  # Internal priming filtering
  IP.seq.I <- paste0(rep("A", times = IA.number), collapse = "")
  results.int <- filter_IP(x = results.int, int.priming.seq = IP.seq.I,
                           left = Ifilter.border.left,
                           right = Ifilter.border.right)
  
  # Statistical analysis ----------------------------------------------------
  # Chi-squre test for APA
  results.int <- test_APA(results.int, clus = "all")
  # Inferring global trends -------------------------------------------------
  results.int <- calc_p_pui_mat(results.int)
  saveRDS(object = results.int, file = "./outs/results_introns.RDS")
  
  # Write final outputs -----------------------------------------------------
  results.int.output <- disply_results(x = results.int, org = org, int = T)
  results.int <- annotate_results(results.int, org = org)
  
    cells_pui <- as.data.frame(results.int@ppui.cells)
    cells_pui$cell <- gsub(pattern = "_proximal_PUI", replacement = "",
                           x = rownames(cells_pui))
    cells_pui$sample <- gsub(pattern = "_.*", replacement = "",
                             x = cells_pui$cell)
    cells_pui$cell <- gsub(pattern = ".*_", replacement = "",
                           x = cells_pui$cell)
    cells_pui <- cells_pui[, c(3, 2, 1)]
    colnames(cells_pui) <- c("Sample", "Cell_BC", "Mean_intronic_PUI")
    write.table(x = cells_pui, file = "./outs/Intronic.Mean.Cell.IPUI.txt",
                quote = F, sep = "\t", col.names = T, row.names = F)
    
  write.table(x = results.int.output, file = "./outs/Intronic.APA.events.txt",
              quote = F, sep = "\t", col.names = T, row.names = F)
  write.table(x = results.int@metadata, file = "./outs/Intronic.peaks.txt",
              quote = F, sep = "\t", col.names = T, row.names = F)
  file.create("./outs/summary.Introns.txt")
  cat("Number of peaks passed fillter:\t", length(results.int@clus.counts),
      "\nNumber of significant (FDR < 5%) APA events:",
      sum(results.int@pvalues[[1]][, 2] < 0.05, na.rm = T),
      "\n", file = "./outs/summary.Introns.txt")
  
  # Plot --------------------------------------------------------------------
  sig.utrs <- as.vector(results.int@pvalues$all[,2] < 0.05)
  sig <- results.int[which(sig.utrs),]
  tidy.ppui <- as.data.frame(sig@ppui.clus)
  colnames(tidy.ppui) <- gsub(x = colnames(tidy.ppui),
                              pattern = "_proximal_PUI", replacement = "")
  tidy.ppui <- tidyr::gather(data = tidy.ppui)
  colnames(tidy.ppui) <- c("Cluster", "value")
  pdf("./outs/Intronic.PUI.ECDF.pdf")
  p <- ggplot2::ggplot(data = tidy.ppui,
                       ggplot2::aes(x = value, color = Cluster))
  p <- p + ggplot2::stat_ecdf(size = 1)
  p <- p + ggplot2::theme_bw()
  p <- p + ggplot2::xlab("Intronic peak usage index")
  p <- p + ggplot2::ylab("Cumulative fraction")
  print(p)
  dev.off()