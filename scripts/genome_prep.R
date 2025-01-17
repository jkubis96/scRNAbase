args <- commandArgs(trailingOnly = TRUE)
set.seed(123)

#Paths and arguments from env
{
  print(args)
  annotation <- args[1]
  output <- args[2]
  source <- args[3]
  three_prime_utr <- as.integer(args[4])
  five_prime_utr <- as.integer(args[5])
  extend <- as.logical(args[6])
  coding_elements <- strsplit(args[7], ",")[[1]]
  shift <- as.integer(args[8])
  optimize <- as.logical(args[9])
  

  
  functions <- file.path(source, 'gtf_tool.R')
  source(functions, local = T)
}

#Configuration file 




GTF <- load_annotation(annotation)



GTF2 <- create_GTF_df(GTF, optimize, shift)

GTF2 <- sort_alias(GTF2)


if (extend && optimize) {

  GTF3 <- add_UTR(GTF2, five_prime_utr, three_prime_utr, coding_elements)
  
} else {
  
  GTF3 <- GTF2
  
}


GTF4 <- create_full_GTF(GTF3)

write.table(GTF4, file.path(output, 'correct_annotation.gtf'), quote = F, sep = '\t', col.names = F, row.names = F)

GTF5 <- create_reduced_GTF(GTF3)


write.table(GTF5, file.path(output, 'reduced_annotation.gtf'), quote = F, sep = '\t', col.names = T, row.names = F)



