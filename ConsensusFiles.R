library(stringr)

args <- commandArgs()

tablename = args[6]
haplohiv_path = args[7]
outpath = args[8]
samples = args[9]

table = read.table(tablename, sep = "\t", header = TRUE)
combi_name = table[grep(str_sub(samples, start=1, end=-7), table$Raw.Filenames), "New.shortname"]
both_samples = table[which(table$New.shortname == combi_name), "Raw.Filenames"]

#write consensus files
names = paste0(haplohiv_path, "check_for_repeats/", str_sub(both_samples, start=1, end=-17), ".fa")
names = str_flatten(names, collapse = " ")
writeLines(names, paste0(outpath, "/consensus/consensusses_",combi_name,".txt"))
