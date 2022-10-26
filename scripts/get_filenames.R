library(stringr)

args <- commandArgs()

combi = args[6]
tablename = args[7]
outpath = args[8]
reads_path = args[9]

table = read.table(tablename, sep = "\t", header = TRUE)
files = table[table$New.shortname == combi, "Filenames"]

#write consensus files
names = paste0(outpath, "/", files, ".fa")
names = str_flatten(names, collapse = " ")
writeLines(names, paste0("names_",combi,".txt"))

#write reads files
readsfiles = paste0(reads_path, "/", files, ".fastq")
readsfiles = str_flatten(readsfiles, collapse = " ")
writeLines(readsfiles, paste0("readfiles_",combi,".txt"))