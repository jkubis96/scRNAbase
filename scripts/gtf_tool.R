load_annotation <- function(path, genetic_elements = NaN) {
  
  set.seed(123)
  
  options(scipen = 999)
  library(readr)
  library(stringr)
  library(dplyr)
  library(doSNOW)
  library(foreach)
  library(doParallel)

  cat('\n\n Data loading... \n\n')


  GTF <- read_delim(path, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, comment = '#', 
                    col_types = cols(
                      X1 = col_character(),
                      X2 = col_character(),
                      X3 = col_character(),
                      X4 = col_integer(),
                      X5 = col_integer(),
                      X6 = col_character(),
                      X7 = col_character(),
                      X8 = col_character(),
                      X9 = col_character()
                    ))
  
  if (!is.na(genetic_elements)) {
    
    
    GTF <- GTF[toupper(GTF$X3) %in% toupper(genetic_elements), ]
    
  }
  
  

return(GTF)

}



sort_alias <- function(input, chromosomes_queue = NaN) {
  
  set.seed(123)
  
  chro_col <- colnames(input)[1]
  start_col <- colnames(input)[4]
  end_col <- colnames(input)[5]
  starnd <- colnames(input)[7]
  
  
  
  if (is.na(chromosomes_queue)) {
    chromosomes_queue = unique(input[[chro_col]])
  }
  
  output <- data.frame()
  
  for (ch in chromosomes_queue) {
    for (s in c('+', '-')) {
      
      cat('\n',paste('CHROMOSOME:', ch , '| STRAND:', s,' sorting...                              '))
      
    
      tmp <- input[(input[[chro_col]] %in% ch & input[[starnd]] %in% s),]
      tmp <- tmp[order(tmp[[start_col]], tmp[[end_col]]),]
      output <- rbind(output, tmp)
      
    }
    
  }
  
  return(output)
}





create_GTF_df <- function(input, optimize = TRUE, shift = 100000) {

    set.seed(123)

    cat('\n\n GTF converting... \n\n')

    df <- input[,1:8]
    
   
    
    if (TRUE %in% unique(grepl('gene_name', input$X9))) {
      
    #GENECODE & https://www.ensembl.org/index.html
    #############################################################
    
      # gene id
      df$gene_id_check <- grepl('gene_id', input$X9)
      
      df$gene_id <- ifelse(
        grepl('gene_id', input$X9),
        gsub(' ', '', gsub('"', '',gsub(".*=", "", gsub(";.*", "", gsub(".*gene_id", "", input$X9))))),
        ""
      )
      
      
      df$gene_id_check <- df$gene_id != ""
      
      
      # gene name
      df$gene_name_check <- grepl('gene_name', input$X9)
      
      df$gene_name <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*=", "", gsub(";.*", "", gsub(".*gene_name", "", input$X9))))),
        ""
      )
      
      
      df$gene_name_check <- df$gene_name != ""
      
      
      
      # transcript name
      df$transcript_name_check <- grepl('transcript_name', input$X9)
      
      df$transcript_name <- ifelse(
        grepl('transcript_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_name", "", input$X9)))))),
        ""
      )
      
      
      df$transcript_name_check <- df$transcript_name != ""
      
      
      
      # transcript id
      df$transcript_id_check <-grepl('transcript_id', input$X9)
      
      df$transcript_id <- ifelse(
        grepl('transcript_id', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*=", "", gsub(".*?\\s", "", gsub(";.*", "", gsub(".*transcript_id", "", input$X9)))))),
        ""
      )
      
      df$transcript_id_check <- df$transcript_id != ""
      
      
    
    if (optimize) {
    
    df$gene_id[df$gene_id == ""] <- df$transcript_id[df$gene_id == ""]
    
    df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE),  ]
    
    df$gene_name[df$gene_name_check == FALSE] <-  df$gene_id[df$gene_name_check == FALSE ]
    

    df$transcript_id[df$transcript_id == ""] <- df$gene_id[df$transcript_id == ""]
    
    df$transcript_name[df$transcript_name_check == FALSE] <- df$gene_name[df$transcript_name_check == FALSE]
    
    #repaire gene names
    gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id', 'transcript_id', 'transcript_name')]
    gene_names <- distinct(gene_names)
   
    
    cat('\n\n Unifying the names of genes... \n\n')
    
    iterations <- length(gene_names$gene_name)
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    df <- df %>% distinct()
    CPU <- detectCores() - 2
    
    cl <- makeCluster(CPU)
    
    
    registerDoParallel(cl)
    registerDoSNOW(cl)
    
    
   
  
    df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
      tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
      tmp$gene_name <- gene_names$gene_name[n]
      tmp
      
    }
    
    close(pb)
    stopCluster(cl)  
    
    
    df <- distinct(df2)
    
    rm(df2)
    
    
    
    
    df$lp <- 1:length(df$X1)
    
    #repair duplicated genes names from different loci
    
    duplicated_names <- unique(df$gene_name[duplicated(df$gene_name)])
    
    
    if (length(duplicated_names) > 0) {
      dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
      dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
      dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
      dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
      
      genes <- unique(duplicated_names)
      global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
      cat('\n\n Duplicated genes repairing...             \n\n ')
      
      for (strand in c('+','-')) {
        tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
        for (gen in genes){
          tmp_ch <- tmp_strand[tmp_strand$gene_name %in% gen, ]
          for (c in unique(tmp_ch$X1)) {
            tmp <- tmp_ch[tmp_ch$X1 %in% c,]
            tmp <- tmp[order(tmp$X4), ]
            group <- 1
            for (i in 1:length(tmp$gene_name)) {
              cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp$gene_name)*100),2), '%           '))
              
              if (length(tmp$gene_name) == 1) {
                
                tmp$group <- group
                global_df <- rbind(global_df, tmp)

              } else if (i  == 1) {
                
                tmp1 <- tmp[i,]
                tmp1$group <- group
                global_df <- rbind(global_df, tmp1)
                
              } else if (i <= length(tmp$gene_name)) {
                tmp1 <- tmp[i-1,]
             
                tmp2 <- tmp[i,]
                if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - shift >  tmp1$X5[1]) {
                  group <- group + 1
                  tmp2$group <- group
                  global_df <- rbind(global_df, tmp2)
                } else {
                  tmp2$group <- group
                  global_df <- rbind(global_df, tmp2)
                }
                
              } 
            }
          } 
        }
      }
      
      
      
      
      
      
      global_df <- global_df %>%
        group_by(renamed) %>%
        filter(mean(group) != 1) %>%
        ungroup()
      
      
      
      
      global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
      
      

      global_df <- global_df %>%
        group_by(gene_name, X1) %>%
        mutate(
          renamed = if (n_distinct(renamed) == 1) {
            substr(renamed, 1, nchar(renamed) - 2)
          } else {
            renamed
          }
        ) %>%
        ungroup()
      
      
      
      
      repaired_df <- data.frame()
      for (g in unique(global_df$renamed)) {
        tmp <- global_df[global_df$renamed %in% g,]
        if (length(unique(tmp$X7)) > 1) {
          tmp <- global_df[global_df$renamed %in% g,]
          tmp$renamed <- paste0(tmp$renamed, '-str.', tmp$X7)
          repaired_df <- rbind(repaired_df, tmp)
        } else {
          repaired_df <- rbind(repaired_df, tmp)
        }
      }
      
      
      
      global_df <- repaired_df
      
      rm(repaired_df)
      
      #rename genes
      
      for (new_name in 1:length(global_df$lp)) {
        cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
        df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
      }
      
      
      
      rm(global_df)
      
    }
    
    }
    
    df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
    
    } else if (TRUE %in% unique(grepl('=gene-', input$X9))) {
      
    #NCBI
    ################################################################
    
     
      
      # gene id
      df$gene_id_check <- grepl('GeneID', input$X9)
      
      df$gene_id <- ifelse(
        grepl('GeneID', input$X9),
        gsub(' ', '', gsub('"', '', gsub(".*?\\s", "", gsub(",.*", "", gsub(".*GeneID:", "", input$X9))))),
        ""
      )
      
      df$gene_id_check <- df$gene_id != ""
      
      
      # gene name 1
      df$gene_name_check <- grepl('=gene-', input$X9)
      
      df$gene_name <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(".*=gene-", "", input$X9)))),
        ""
      )
      
      
      df$gene_name_check <- df$gene_name != ""
      
      
      # gene name 1
      df$gene_name_check2 <- grepl(';gene=', input$X9)
      
      df$gene_name2 <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(".*;gene=", "", input$X9)))),
        ""
      )
      
      
      df$gene_name_check2 <- df$gene_name2 != ""
      
      # transcript name
      
      df$transcript_name <- NA
      
      
      # transcript id
      df$transcript_id_check <-grepl('Genbank:', input$X9)
      
      df$transcript_id <- ifelse(
        grepl('transcript_id', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(",.*", "", gsub(".*Genbank:", "", input$X9))))),
        ""
      )
      
      df$transcript_id_check <- df$transcript_id != ""
      
      
      if (optimize) {
      #
      
      
      df$gene_id[df$gene_id == ""] <- df$transcript_id[df$gene_id == ""]
      
      df$gene_name[df$gene_name_check == FALSE & df$gene_name_check2 == TRUE] <- df$gene_name2[df$gene_name_check2 == TRUE & df$gene_name_check == FALSE]
      
      df$gene_name_check[df$gene_name_check2 == TRUE] <- TRUE
      
      df <- df[!(df$gene_name_check == FALSE & df$gene_id_check == FALSE & df$transcript_id_check == FALSE),  ]
      
      df$gene_name[df$gene_name_check == FALSE] <-  df$gene_id[df$gene_name_check == FALSE ]
      
      df$transcript_id[df$transcript_id == ""] <- df$gene_id[df$transcript_id == ""]
      
      df$transcript_name <- df$transcript_id
      
      
      #repaire gene names
      gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id', 'transcript_id', 'transcript_name')]
      gene_names <- distinct(gene_names)
      
      
      cat('\n\n Unifying the names of genes... \n\n')
      
      iterations <- length(gene_names$gene_name)
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      df <- df %>% distinct()
      CPU <- detectCores() - 2
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
      
      
      
      
      df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
        tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
        tmp$gene_name <- gene_names$gene_name[n]
        tmp
        
      }
      
      close(pb)
      stopCluster(cl)  
      
      
      df <- distinct(df2)
      
      rm(df2)
      
      
      
      
      df$lp <- 1:length(df$X1)
      
      #repair duplicated genes names from different loci
      
      duplicated_names <- unique(df$gene_name[duplicated(df$gene_name)])
      
      
      if (length(duplicated_names) > 0) {
        dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
        dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
        dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
        dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
        
        genes <- unique(duplicated_names)
        global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
        cat('\n\n Duplicated genes repairing...             \n\n ')
        
        for (strand in c('+','-')) {
          tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
          for (gen in genes){
            tmp_ch <- tmp_strand[tmp_strand$gene_name %in% gen, ]
            for (c in unique(tmp_ch$X1)) {
              tmp <- tmp_ch[tmp_ch$X1 %in% c,]
              tmp <- tmp[order(tmp$X4), ]
              group <- 1
              for (i in 1:length(tmp$gene_name)) {
                cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp$gene_name)*100),2), '%           '))
                
                if (length(tmp$gene_name) == 1) {
                  
                  tmp$group <- group
                  global_df <- rbind(global_df, tmp)
                  
                } else if (i  == 1) {
                  
                  tmp1 <- tmp[i,]
                  tmp1$group <- group
                  global_df <- rbind(global_df, tmp1)
                  
                } else if (i <= length(tmp$gene_name)) {
                  tmp1 <- tmp[i-1,]
                  
                  tmp2 <- tmp[i,]
                  if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - shift >  tmp1$X5[1]) {
                    group <- group + 1
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  } else {
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  }
                  
                } 
              }
            } 
          }
        }
        
        
        
        
        
        
        global_df <- global_df %>%
          group_by(renamed) %>%
          filter(mean(group) != 1) %>%
          ungroup()
        
        
        
        
        global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
        
        
        
        global_df <- global_df %>%
          group_by(gene_name, X1) %>%
          mutate(
            renamed = if (n_distinct(renamed) == 1) {
              substr(renamed, 1, nchar(renamed) - 2)
            } else {
              renamed
            }
          ) %>%
          ungroup()
        
        
        
        
        repaired_df <- data.frame()
        for (g in unique(global_df$renamed)) {
          tmp <- global_df[global_df$renamed %in% g,]
          if (length(unique(tmp$X7)) > 1) {
            tmp <- global_df[global_df$renamed %in% g,]
            tmp$renamed <- paste0(tmp$renamed, '-str.', tmp$X7)
            repaired_df <- rbind(repaired_df, tmp)
          } else {
            repaired_df <- rbind(repaired_df, tmp)
          }
        }
        
        
        
        global_df <- repaired_df
        
        rm(repaired_df)
        
        #rename genes
        
        for (new_name in 1:length(global_df$lp)) {
          cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
          df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
        }
        
        
        
        rm(global_df)
        
      }
      
      }
      
      df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
      
  
    } else if (TRUE %in% unique(grepl('gene:', input$X9)) & TRUE %in% unique(grepl('transcript:', input$X9))) {
      
      #CUSTOM
      ################################################################
      
    

      # gene name
      df$gene_name_check <- grepl('gene:', input$X9)
      
      df$gene_name <- ifelse(
        grepl('gene_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(".*gene:", "", input$X9)))),
        ""
      )
      
      
      df$gene_name_check <- df$gene_name != ""
      
      # transcript name
      df$transcript_name_check <- grepl('transcript:', input$X9)
      
      df$transcript_name <- ifelse(
        grepl('transcript_name', input$X9),
        gsub(' ', '', gsub('"', '', gsub(";.*", "", gsub(",.*", "", gsub(".*transcript:", "", input$X9))))),
        ""
      )
      

      
      # transcript id
      # gene id
      
      
      df$transcript_id <- df$transcript_name
      
      df$gene_id <- df$gene_name
      
      
      if (optimize) {
        
      
      df <- df[!(df$gene_name_check == FALSE & df$transcript_name_check == FALSE),  ]
      
      df$gene_name[df$gene_name_check == FALSE] <-  df$transcript_name[df$gene_name_check == FALSE ]
      df$transcript_name[df$transcript_name_check == FALSE] <-  df$gene_name[df$transcript_name_check == FALSE ]
      
      #
      
   
      
      
        
      #repaire gene names
      gene_names <- df[,colnames(df) %in% c('gene_name', 'gene_id', 'transcript_id', 'transcript_name')]
      gene_names <- distinct(gene_names)
      
      
      cat('\n\n Unifying the names of genes... \n\n')
      
      iterations <- length(gene_names$gene_name)
      pb <- txtProgressBar(max = iterations, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      
      df <- df %>% distinct()
      CPU <- detectCores() - 2
      
      cl <- makeCluster(CPU)
      
      
      registerDoParallel(cl)
      registerDoSNOW(cl)
      
      
      
      
      df2 <- foreach(n = 1:length(gene_names$gene_name), .options.snow = opts, .combine=rbind) %dopar% {
        tmp <- df[df$gene_id %in% gene_names$gene_id[n],] 
        tmp$gene_name <- gene_names$gene_name[n]
        tmp
        
      }
      
      close(pb)
      stopCluster(cl)  
      
      
      df <- distinct(df2)
      
      rm(df2)
      
      
      
      
      df$lp <- 1:length(df$X1)
      
      #repair duplicated genes names from different loci
      
      duplicated_names <- unique(df$gene_name[duplicated(df$gene_name)])
      
      
      if (length(duplicated_names) > 0) {
        dedup_df <- df[, colnames(df) %in% c('X1','X7', 'X3', 'X4', 'X5', 'gene_name', 'lp')]
        dedup_df <- dedup_df[dedup_df$gene_name %in% duplicated_names,]
        dedup_df$renamed <- gsub(' ', '', gsub('"', '', dedup_df$gene_name))
        dedup_df$renamed <- paste0(dedup_df$renamed, '.',  dedup_df$X1)
        
        genes <- unique(duplicated_names)
        global_df <- data.frame(matrix(ncol = length(colnames(dedup_df)), nrow = 0))
        cat('\n\n Duplicated genes repairing...             \n\n ')
        
        for (strand in c('+','-')) {
          tmp_strand <- dedup_df[dedup_df$X7 %in% strand,]
          for (gen in genes){
            tmp_ch <- tmp_strand[tmp_strand$gene_name %in% gen, ]
            for (c in unique(tmp_ch$X1)) {
              tmp <- tmp_ch[tmp_ch$X1 %in% c,]
              tmp <- tmp[order(tmp$X4), ]
              group <- 1
              for (i in 1:length(tmp$gene_name)) {
                cat('\r',paste('GENE:',gen ,'| STRAND:', strand , '| PROGRESS:', round((i/length(tmp$gene_name)*100),2), '%           '))
                
                if (length(tmp$gene_name) == 1) {
                  
                  tmp$group <- group
                  global_df <- rbind(global_df, tmp)
                  
                } else if (i  == 1) {
                  
                  tmp1 <- tmp[i,]
                  tmp1$group <- group
                  global_df <- rbind(global_df, tmp1)
                  
                } else if (i <= length(tmp$gene_name)) {
                  tmp1 <- tmp[i-1,]
                  
                  tmp2 <- tmp[i,]
                  if (tmp2$X1[1] == tmp1$X1[1] & tmp2$X4[1] - shift >  tmp1$X5[1]) {
                    group <- group + 1
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  } else {
                    tmp2$group <- group
                    global_df <- rbind(global_df, tmp2)
                  }
                  
                } 
              }
            } 
          }
        }
        
        
        
        
        
        
        global_df <- global_df %>%
          group_by(renamed) %>%
          filter(mean(group) != 1) %>%
          ungroup()
        
        
        
        
        global_df$renamed <- paste0(global_df$renamed, '.', global_df$group)
        
        
        
        global_df <- global_df %>%
          group_by(gene_name, X1) %>%
          mutate(
            renamed = if (n_distinct(renamed) == 1) {
              substr(renamed, 1, nchar(renamed) - 2)
            } else {
              renamed
            }
          ) %>%
          ungroup()
        
        
        
        
        repaired_df <- data.frame()
        for (g in unique(global_df$renamed)) {
          tmp <- global_df[global_df$renamed %in% g,]
          if (length(unique(tmp$X7)) > 1) {
            tmp <- global_df[global_df$renamed %in% g,]
            tmp$renamed <- paste0(tmp$renamed, '-str.', tmp$X7)
            repaired_df <- rbind(repaired_df, tmp)
          } else {
            repaired_df <- rbind(repaired_df, tmp)
          }
        }
        
        
        
        global_df <- repaired_df
        
        rm(repaired_df)
        
        #rename genes
        
        for (new_name in 1:length(global_df$lp)) {
          cat('\r',paste('Gene renaming -> PROGRESS:', round((new_name/length(global_df$lp)*100),2), '%                            '))
          df$gene_name[df$lp == global_df$lp[new_name]] <- global_df$renamed[global_df$lp == global_df$lp[new_name]]
        }
        
        
        
        rm(global_df)
        
      }
      
      }
      
      df <- df[,!colnames(df) %in% c('lp', 'gene_id_check', 'gene_name_check', 'transcript_name_check', 'transcript_id_check')]
      
      
    }
    
    
    
    colnames(df) <- c('chr','source','annotationType','start','end','score','strand','phase','gene_id','gene_name','transcript_name','transcript_id')

return(df)

}




add_UTR <- function(input, five_prime_utr_length = 400, three_prime_utr_length = 800, genetic_elements = c("EXON", "CDS", 'TRANSCRIPT', 'MRNA')) {
  
  set.seed(123)

  cat('\n\n UTRs sequence extending...             \n\n')
  
  tmp_final <- data.frame()
  final <- input
  
  # input$sort_val <- 1:length(input$chr)
  chromosomes <- unique(input$chr)
  
  CDS <- final[toupper(input$annotationType) %in% toupper(genetic_elements), ]
  
  
  for (chr in chromosomes) {
  
  tmp_primary <- CDS[CDS$chr %in% chr,]
  # tmp_primary <- tmp_primary[order(tmp_primary$start),]
  # tmp_primary$lp <- 1:length(tmp_primary$chr)
  
  
  for (sen in c('+','-')) {
  tmp <- tmp_primary[tmp_primary$strand %in% sen,]
  gen_list <- unique(tmp$gene_name)

  
  if (length(gen_list) > 0) {
    
  
  for (i in 1:length(gen_list)) {
          # debug
          # if (i == 12) {break}
  
          # tmp variables
    
          five_prime_utr_length_cor_minus = NaN
          five_prime_utr_length_cor_plus = NaN
          three_prime_utr_length_cor_minus = NaN
          three_prime_utr_length_cor_plus = NaN
            
          cat('\r',paste('LOC:',chr, '| STRAND:', sen , '| PROGRESS:', round((i/length(gen_list)*100),2), '%           '))
          
          if (i == 1) {
            
            # first element
            
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            
            tmp1_min <- tmp1[tmp1$start == min(tmp1$start),]
            tmp1_min <- tmp1_min[1,]
            tmp1_max <- tmp1[tmp1$end == max(tmp1$end),]
            tmp1_max <- tmp1_max[1,]
            
            
            tmp2 <- tmp[tmp$gene_name %in% gen_list[i+1],] 
            
            
            if (min(tmp2$start) > (max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)+2)) {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length,0)
              } else {
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
              }
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/2) {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.5,0)
              } else {
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.5,0)
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
              }
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/3) {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
                three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.3,0)
              } else {
                five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.3,0)
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                
              }
              
            } else if (min(tmp2$start) > max(tmp1$end)){
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- five_prime_utr_length
                three_prime_utr_length_cor_plus <- round((min(tmp2$start) - max(tmp1$end) -2) / 2,0)
                if (three_prime_utr_length_cor_plus < 0) {
                  three_prime_utr_length_cor_plus <- 0
                }
              } else {
                five_prime_utr_length_cor_minus <- round((min(tmp2$start) - max(tmp1$end) -2) / 2, 0)
                if (five_prime_utr_length_cor_minus < 0) {
                  five_prime_utr_length_cor_minus <- 0
                }
                three_prime_utr_length_cor_minus <- three_prime_utr_length
                
              }
              
              
            } else {
              
              if (sen == '+') {
                five_prime_utr_length_cor_plus <- five_prime_utr_length
               
              } else {
                
                three_prime_utr_length_cor_minus <- three_prime_utr_length
                
              }
              
              
            }
            
            
            
            
            tmp0 <- tmp1
            tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
            tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
            tmp0 <- tmp0[1,]
            
            tmp0_UTR5 <- tmp0
            tmp0_UTR5$source <- 'JBIO-predicted'
            tmp0_UTR5$annotationType <- 'five_prime_UTR'
            tmp0_UTR5$score <- '.'
            tmp0_UTR5$phase <- '.'
            
            
            tmp0_UTR3 <- tmp0
            tmp0_UTR3$source <- 'JBIO-predicted'
            tmp0_UTR3$annotationType <- 'three_prime_UTR'
            tmp0_UTR3$score <- '.'
            tmp0_UTR3$phase <- '.'
            
            tmp0_transcript <- tmp0
            tmp0_transcript$source <- 'JBIO-predicted'
            tmp0_transcript$annotationType <- 'transcript'
            tmp0_transcript$score <- '.'
            tmp0_transcript$phase <- '.'
            
            if (sen == '+') {
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_plus)) {
                
                utr5_start <- as.integer(tmp1_min$start - five_prime_utr_length_cor_plus)
                if (utr5_start < 0) {
                  utr5_start <- 0
                }
                
                tmp0_UTR5$start <- utr5_start
                tmp0_UTR5$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR3$end <- as.integer(tmp1_max$end + three_prime_utr_length_cor_plus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #TRANSCRIPT
              
              if (!is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
              tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
              tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
              
              tmp_final <- rbind(tmp_final, tmp0_transcript)
              tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(three_prime_utr_length_cor_plus) & is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
              
            } else if (sen == '-') {
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_minus)) {
                
                utr3_start <- as.integer(tmp1_min$start - three_prime_utr_length_cor_minus)
                if (utr3_start < 0) {
                  utr3_start <- 0
                }
                
                tmp0_UTR3$start <- utr3_start
                tmp0_UTR3$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_minus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR5$end <- as.integer(tmp1_max$end + five_prime_utr_length_cor_minus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
            }
            
            
          } else if (i < length(gen_list)) {
            
            # middle element
            
            tmp_p <- tmp[tmp$gene_name %in% gen_list[i-1],]
            
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            
            
            tmp1_min <- tmp1[tmp1$start == min(tmp1$start),]
            tmp1_min <- tmp1_min[1,]
            tmp1_max <- tmp1[tmp1$end == max(tmp1$end),]
            tmp1_max <- tmp1_max[1,]
            
            
            tmp2 <- tmp[tmp$gene_name %in% gen_list[i+1],] 
            
            
            if (min(tmp2$start) > (max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)+2)) {
              
              three_prime_utr_length_cor_plus <- three_prime_utr_length
      
              five_prime_utr_length_cor_minus <- five_prime_utr_length
              
              ################################################################################################
              if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
                
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                  
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/2) {
              
              three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.5,0)
        
              five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.5,0)
              
              ################################################################################################
              if (max(tmp_p$end) > (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) > (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                  
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
              
              
            } else if (min(tmp2$start) > max(tmp1$end) + (five_prime_utr_length + three_prime_utr_length)/3) {
              
              three_prime_utr_length_cor_plus <- round(three_prime_utr_length*0.3,0)
              
              five_prime_utr_length_cor_minus <- round(five_prime_utr_length*0.3,0)
              
              ################################################################################################
              if (max(tmp_p$end) > (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) + (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
              
            } else if (min(tmp2$start) > max(tmp1$end)) {
              
            
              
              three_prime_utr_length_cor_plus <- round((min(tmp2$start) - max(tmp1$end) - 2)/2,0)
                if (three_prime_utr_length_cor_plus < 0) {
                  three_prime_utr_length_cor_plus <- 0
                }
              
              five_prime_utr_length_cor_minus <- round((min(tmp2$start) - max(tmp1$end) - 2)/2,0)
                if (five_prime_utr_length_cor_minus < 0) {
                  five_prime_utr_length_cor_minus <- 0
                }

              ################################################################################################
              if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              #################################################################################################
              
            } else if (max(tmp_p$end) < min(tmp1$start)) {
              
              
              if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)+2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length )/2) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length )/2) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
                #
              } else if (max(tmp_p$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              } else if (max(tmp_p$end) < min(tmp1$start)) {
                
                if (sen == '+') {
                  five_prime_utr_length_cor_plus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (five_prime_utr_length_cor_plus < 0) {
                    five_prime_utr_length_cor_plus <- 0
                  }
                } else {
                  three_prime_utr_length_cor_minus <- min(tmp1$start) - max(tmp_p$end) -2
                  if (three_prime_utr_length_cor_minus < 0) {
                    three_prime_utr_length_cor_minus <- 0
                  }
                  
                }
                
              }
              
              
            }
            
            
            tmp0 <- tmp1
            tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
            tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
            tmp0 <- tmp0[1,]
            
            tmp0_UTR5 <- tmp0
            tmp0_UTR5$source <- 'JBIO-predicted'
            tmp0_UTR5$annotationType <- 'five_prime_UTR'
            tmp0_UTR5$score <- '.'
            tmp0_UTR5$phase <- '.'
            
            
            tmp0_UTR3 <- tmp0
            tmp0_UTR3$source <- 'JBIO-predicted'
            tmp0_UTR3$annotationType <- 'three_prime_UTR'
            tmp0_UTR3$score <- '.'
            tmp0_UTR3$phase <- '.'
            
            tmp0_transcript <- tmp0
            tmp0_transcript$source <- 'JBIO-predicted'
            tmp0_transcript$annotationType <- 'transcript'
            tmp0_transcript$score <- '.'
            tmp0_transcript$phase <- '.'
            
            if (sen == '+') {
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_min$start - five_prime_utr_length_cor_plus)
                tmp0_UTR5$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR3$end <- as.integer(tmp1_max$end + three_prime_utr_length_cor_plus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #TRANSCRIPT
              if (!is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (is.na(three_prime_utr_length_cor_plus) & !is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(three_prime_utr_length_cor_plus) & is.na(five_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
              
              
            } else if (sen == '-') {
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_min$start - three_prime_utr_length_cor_minus)
                tmp0_UTR3$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_minus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR5$end <- as.integer(tmp1_max$end + five_prime_utr_length_cor_minus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
                
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
            }
            
          } else if (i == length(gen_list)) {
            
            # final element
            
            tmp1 <- tmp[tmp$gene_name %in% gen_list[i],]
            
            tmp1_min <- tmp1[tmp1$start == min(tmp1$start),]
            tmp1_min <- tmp1_min[1,]
            tmp1_max <- tmp1[tmp1$end == max(tmp1$end),]
            tmp1_max <- tmp1_max[1,]
            
            tmp2 <- tmp[tmp$gene_name %in% gen_list[i-1],] 
            
            
            # End UTR
            
            if (max(tmp2$end) < max(tmp1$end)) {
              
              three_prime_utr_length_cor_plus <- round(three_prime_utr_length,0)
              five_prime_utr_length_cor_minus <- round(five_prime_utr_length,0)
              
            }
              
            
            # tmp2 -> tmp1 +
            if (max(tmp2$end) < (min(tmp1$start) - (five_prime_utr_length )+2) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (five_prime_utr_length)/2) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.5,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (five_prime_utr_length)/3) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round(five_prime_utr_length*0.3,0)
              
            } else if (max(tmp2$end) < min(tmp1$start) & sen == '+') {
              
                five_prime_utr_length_cor_plus <- round((min(tmp1$start) - max(tmp2$end) -2), 0)
                if (five_prime_utr_length_cor_plus < 0) {
                  five_prime_utr_length_cor_plus <- 0
                }
              
            }
            
            
            # tmp2 -> tmp1 -
            if (max(tmp2$end) < (min(tmp1$start) - (three_prime_utr_length)+2) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (three_prime_utr_length)/2) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.5,0)
              
            } else if (max(tmp2$end) < (min(tmp1$start) - (three_prime_utr_length)/3) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round(three_prime_utr_length*0.3,0)
              
            } else if (max(tmp2$end) < min(tmp1$start) & sen == '-') {
              
                three_prime_utr_length_cor_minus <- round((min(tmp1$start) - max(tmp2$end) -2) / 2, 0)
                if (three_prime_utr_length_cor_minus < 0) {
                  three_prime_utr_length_cor_minus <- 0
               
                }
              
            }
            
            
            
            tmp0 <- tmp1
            tmp0 <- tmp0[tmp0$start == min(tmp0$start),]
            tmp0 <- tmp0[tmp0$end == max(tmp0$end),]
            tmp0 <- tmp0[1,]
            
            tmp0_UTR5 <- tmp0
            tmp0_UTR5$source <- 'JBIO-predicted'
            tmp0_UTR5$annotationType <- 'five_prime_UTR'
            tmp0_UTR5$score <- '.'
            tmp0_UTR5$phase <- '.'
            
            
            tmp0_UTR3 <- tmp0
            tmp0_UTR3$source <- 'JBIO-predicted'
            tmp0_UTR3$annotationType <- 'three_prime_UTR'
            tmp0_UTR3$score <- '.'
            tmp0_UTR3$phase <- '.'
            
            tmp0_transcript <- tmp0
            tmp0_transcript$source <- 'JBIO-predicted'
            tmp0_transcript$annotationType <- 'transcript'
            tmp0_transcript$score <- '.'
            tmp0_transcript$phase <- '.'
            
            if (sen == '+') {
              
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_plus)) {
                
                utr5_start <- as.integer(tmp1_min$start - five_prime_utr_length_cor_plus)
                if (utr5_start < 0) {
                  utr5_start <- 0
                }
                
                tmp0_UTR5$start <- utr5_start
                tmp0_UTR5$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_UTR3$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR3$end <- as.integer(tmp1_max$end + three_prime_utr_length_cor_plus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
              
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_plus) & !is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(five_prime_utr_length_cor_plus) & !is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR3$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_plus) & is.na(three_prime_utr_length_cor_plus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR5$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
              
              
              
            } else if (sen == '-') {
              
              #UTR3
              if (!is.na(three_prime_utr_length_cor_minus)) {
                
                utr3_start <- as.integer(tmp1_min$start - three_prime_utr_length_cor_minus)
                if (utr3_start < 0) {
                  utr3_start <- 0
                }
              
                tmp0_UTR3$start <- utr3_start
                tmp0_UTR3$end <- as.integer(tmp1_min$start - 1)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR3)

              }
                
              
              #UTR5
              if (!is.na(five_prime_utr_length_cor_minus)) {
                
                tmp0_UTR5$start <- as.integer(tmp1_max$end + 1)
                tmp0_UTR5$end <- as.integer(tmp1_max$end + five_prime_utr_length_cor_minus)
                
                tmp_final <- rbind(tmp_final, tmp0_UTR5)

              }
              
              
              #TRANSCRIPT
              if (!is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
              
              } else if (is.na(five_prime_utr_length_cor_minus) & !is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp0_UTR3$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else if (!is.na(five_prime_utr_length_cor_minus) & is.na(three_prime_utr_length_cor_minus)) {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp0_UTR5$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              } else {
                
                tmp0_transcript$start <- as.integer(tmp1_min$start)
                tmp0_transcript$end <- as.integer(tmp1_max$end)
                
                tmp_final <- rbind(tmp_final, tmp0_transcript)
                tmp <- rbind(tmp, tmp0_transcript)
                
              }
                
            }
            
            
          } 
          

        } 
      }
    }
  }

  
  final <- rbind(final, tmp_final)
  
  output <- sort_alias(final)
  
 
  
  return(output)
}







refflat_create <- function(input, geneName = 'gene_name', name = 'gene_id') {
  
  set.seed(123)

  iterations <- length(unique(input$chr))
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  CPU <- detectCores() - 2
  
  cl <- makeCluster(CPU)

  
  
  registerDoParallel(cl)
  registerDoSNOW(cl)
  
  
  
  
  input <- input[,c('chr', 'start', 'end' ,'gene_name', 'gene_id', 'transcript_name', 'transcript_id', 'strand', 'annotationType')]
  input$annotationType[toupper(input$annotationType) %in% c('MRNA', 'GENE', 'TRANSCRIPT')] <- 'TRANSCRIPT'
  input <- distinct(input)
  
  input$sort_val <- 1:length(input$chr)
  input$diff <- input$end - input$start
  

  chromosomes <- unique(input$chr)

  
  df <- data.frame(matrix(ncol = 11, nrow = 0))
  colnames(df) <- c('geneName', 'name', 'chrom', 'strand', 'txStart','txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds')
  
  cat('\n\n Refflat creating ... \n\n')
  results <- foreach(chr = chromosomes, .packages = c('dplyr'), .options.snow = opts, .combine=rbind) %dopar% {
    
    tmp_primary <- input[input$chr %in% chr,]
    tmp_primary <- tmp_primary[order(tmp_primary$start),]
    tmp_primary$lp <- 1:length(tmp_primary$chr)
    
    
    for (sen in c('+','-')) {
      tmp <- tmp_primary[tmp_primary$strand %in% sen,]
      gen_list <- unique(tmp$gene_name)
      
      if (length(gen_list) == 1) {break}
      
      for (i in 1:length(gen_list)) {
        

          check_vector <- toupper(tmp$annotationType[tmp$gene_name %in% gen_list[i]])
          if ('EXON' %in% check_vector & 'TRANSCRIPT' %in% check_vector) {
          tmp_tx <- tmp[tmp$gene_name %in% gen_list[i],]
          tmp_tx <- tmp_tx[toupper(tmp_tx$annotationType) == 'TRANSCRIPT',]
          tmp_dif <- tmp_tx[tmp_tx$diff == max(tmp_tx$diff),]
          min_tr <- min(tmp_dif$start)
          max_tr <- max(tmp_dif$end)
          tmp_tx <- tmp_tx[tmp_tx$diff == max(tmp_tx$diff) | tmp_tx$start < min_tr | tmp_tx$end > max_tr ,]
          tmp_ex <- tmp[tmp$gene_name %in% gen_list[i] & toupper(tmp$annotationType) == 'EXON',]

          for (trancript in 1:nrow(tmp_tx)) {
            tmp_tx_tmp <- tmp_tx[trancript,]
            tmp_ex_tmp <- tmp_ex[tmp_ex$start >= tmp_tx_tmp$start & tmp_ex$end <= tmp_tx_tmp$end,]
            
            decission <- TRUE
            trn <- 0
            while (decission) {
              trn <- trn + 1
                dec <- c()
                start <- c()
                end <- c()
                before <- 0
                for(ex in 1:nrow(tmp_ex_tmp)) {
                  if (before <= tmp_ex_tmp$start[ex]) {
                    start <- c(start, as.integer(tmp_ex_tmp$start[ex]))
                    end <- c(end, as.integer(tmp_ex_tmp$end[ex]))
                    before  <-  as.integer(tmp_ex_tmp$end[ex])
                    
                   
      
                  } else if (before > as.integer(tmp_ex_tmp$start[ex])) {
                    
                    dec <- c(dec, ex-1)
                   
                  }
                }
                
                if (length(dec) > 0) {
                tmp_ex_tmp <- tmp_ex_tmp[-dec[1], ]
                decission <- TRUE
                
                } else if (length(dec) == 0) {
                  decission <- FALSE
                  
                }
                
                
                df[nrow(df) + 1,] <- c(as.character(gsub('"', '', tmp_tx_tmp[[geneName]][1])), as.character(gsub('"', '', tmp_ex_tmp[[name]][1])), as.character(tmp_ex_tmp$chr[1]), as.character(tmp_ex_tmp$strand[1]), as.integer(tmp_tx_tmp$start[1]), as.integer(tmp_tx_tmp$end[1]), as.integer(min(start)), as.integer(max(end)), as.integer(length(start)), sub('"', '', paste0(start, collapse = ',')), sub('"', '', paste0(end, collapse = ',')))
                df <- distinct(df)
                
                
                } 
                  
                 
                
                }
              
            
                }
          
          
              }
            }
                results <- df
        }
  
      
      close(pb)
      stopCluster(cl)  
      return(results)
} 





###############

create_full_GTF <- function(input) {
  
  set.seed(123)

  output <- input[,1:8]
  gene_id <- paste0('gene_id ',gsub(' ', '',input$gene_id), ';')
  gene_name <-  paste0('gene_name ', gsub(' ','',input$gene_name), ';')
  transcript_name <- paste0('transcript_name ', gsub(' ', '',input$transcript_name), ';')
  transcript_id <- paste0('transcript_id ', gsub(' ', '',input$transcript_id), ';')
  output$combine <- paste0(gene_id, gene_name, transcript_name, transcript_id)
  
  
  return(output)
  
}




create_reduced_GTF <- function(input) {
  
  set.seed(123)

  output <- input %>%
    dplyr::select(chr,start,end,strand,transcript_name,transcript_id, gene_name, gene_id)
  output$annotationType <- input$annotationType
  output$transcriptType <- input$source
  
  
  return(output)  
  
}

