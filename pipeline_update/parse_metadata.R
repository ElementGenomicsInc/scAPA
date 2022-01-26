df = read.delim("/data/GSE114530_metadata.tsv")

reformat = function(string){
  return(strsplit(string, "_")[[1]][1])
}

reformat2 = function(string){
  return(gsub("[/ ]", "_", string))
}

df = df[which(df$orig.ident == "week16"),]

df2 = data.frame("source_classification"=unlist(lapply(df$source_classification, reformat2)), 
                 row.names = lapply(rownames(df), reformat))
write.table(df2, "/data/GSE114530/scAPA_input/clusters_SRR7171583.txt", sep="\t", col.names=FALSE, 
            quote=FALSE)
