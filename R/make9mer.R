#' make9mer
#'
#' 9mer generator from protein sequence
#' @details needs a protein length at least 18 amino acid
#' @param protseq Sequence of the protein. Data.frame with: codon number in codon column, and amino acid for each codon in aminoacid column
#' 
#' @return list9mers. List of length equal of protein length. Includes data frames with 9-mers.
#' 
#' @author Nathan Lemonnier \email{nathanael.lemonnier@@gmail.com}
#' 
#' @examples
#' protseq <- data.frame(  codon = seq(1:20)
#'                       , aminoacid = c("M","E","E","P","Q","S","D","P","S","V","E","P","P","L","S","Q","E","T","F","S","D")
#'                       )
#' make9mer(protseq = protseq)

make9mer <- function(protseq){

cat("...creating directories...\n")
#create directories
dir.create(file.path("./", "./files"))
  
  
  AAseq   <- as.data.frame(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"))
  colnames(AAseq)[1] <- c("V1")
  
  protseq <- as.character(protseq$aminoacid)
  codon <- length(protseq)
  AAn <- nrow(AAseq)
  
  list9mers <- vector("list", codon)
  
  # 9mers CODE
  #----
  
  for(i in 1:codon){
    
    #codon=1  
    if(i==1){
      table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
      colnames(table9mers) <- AAseq$V1
      
      for (a in i){ 
        for (b in 1:ncol(table9mers)){
          table9mers[i,b] <- paste(colnames(table9mers)[b]
                                   , protseq[2]
                                   , protseq[3]
                                   , protseq[4]
                                   , protseq[5]
                                   , protseq[6]
                                   , protseq[7]
                                   , protseq[8]
                                   , protseq[9]
                                   ,sep="")}}
      list9mers[[i]] <- table9mers
    }else{
      
      #codon=2  
      if(i==2){
        table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
        colnames(table9mers) <- AAseq$V1
        for (a in 1:i){       
          for (b in 1:ncol(table9mers)){
            if(a==i-1){
              table9mers[a,b] <- paste(colnames(table9mers)[b]
                                       , protseq[i+1]
                                       , protseq[i+2]
                                       , protseq[i+3]
                                       , protseq[i+4]
                                       , protseq[i+5]
                                       , protseq[i+6]
                                       , protseq[i+7]
                                       , protseq[i+8]
                                       ,sep="")
            }else{if(a==i){
              table9mers[a,b] <- paste(protseq[i-1]
                                       , colnames(table9mers)[b]
                                       , protseq[i+1]
                                       , protseq[i+2]
                                       , protseq[i+3]
                                       , protseq[i+4]
                                       , protseq[i+5]
                                       , protseq[i+6]
                                       , protseq[i+7]
                                       ,sep="")}}}}
        list9mers[[i]] <- table9mers  
      }else{
        
        #codon=3           
        if(i==3){
          table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
          colnames(table9mers) <- AAseq$V1
          for (a in 1:i){       
            for (b in 1:ncol(table9mers)){
              if(a==i-2){
                table9mers[a,b] <- paste(colnames(table9mers)[b]
                                         , protseq[i+1]
                                         , protseq[i+2]
                                         , protseq[i+3]
                                         , protseq[i+4]
                                         , protseq[i+5]
                                         , protseq[i+6]
                                         , protseq[i+7]
                                         , protseq[i+8]
                                         ,sep="")
              }else{if(a==i-1){
                table9mers[a,b] <- paste(protseq[i-1]
                                         , colnames(table9mers)[b]
                                         , protseq[i+1]
                                         , protseq[i+2]
                                         , protseq[i+3]
                                         , protseq[i+4]
                                         , protseq[i+5]
                                         , protseq[i+6]
                                         , protseq[i+7]
                                         ,sep="")
              }else{if(a==i){
                table9mers[a,b] <- paste(protseq[i-2]
                                         , protseq[i-1]
                                         , colnames(table9mers)[b]
                                         , protseq[i+1]
                                         , protseq[i+2]
                                         , protseq[i+3]
                                         , protseq[i+4]
                                         , protseq[i+5]
                                         , protseq[i+6]
                                         ,sep="")}}}}}
          list9mers[[i]] <- table9mers
        }else{
          
          #codon=4           
          if(i==4){
            table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
            colnames(table9mers) <- AAseq$V1
            for (a in 1:i){       
              for (b in 1:ncol(table9mers)){
                if(a==i-3){
                  table9mers[a,b] <- paste(colnames(table9mers)[b]
                                           , protseq[i+1]
                                           , protseq[i+2]
                                           , protseq[i+3]
                                           , protseq[i+4]
                                           , protseq[i+5]
                                           , protseq[i+6]
                                           , protseq[i+7]
                                           , protseq[i+8]
                                           ,sep="")
                }else{if(a==i-2){
                  table9mers[a,b] <- paste(protseq[i-1]
                                           , colnames(table9mers)[b]
                                           , protseq[i+1]
                                           , protseq[i+2]
                                           , protseq[i+3]
                                           , protseq[i+4]
                                           , protseq[i+5]
                                           , protseq[i+6]
                                           , protseq[i+7]
                                           ,sep="")
                }else{if(a==i-1){
                  table9mers[a,b] <- paste(protseq[i-2]
                                           , protseq[i-1]
                                           , colnames(table9mers)[b]
                                           , protseq[i+1]
                                           , protseq[i+2]
                                           , protseq[i+3]
                                           , protseq[i+4]
                                           , protseq[i+5]
                                           , protseq[i+6]
                                           ,sep="")
                }else{if(a==i){
                  table9mers[a,b] <- paste(protseq[i-3]
                                           , protseq[i-2]
                                           , protseq[i-1]
                                           , colnames(table9mers)[b]
                                           , protseq[i+1]
                                           , protseq[i+2]
                                           , protseq[i+3]
                                           , protseq[i+4]
                                           , protseq[i+5]
                                           ,sep="")}}}}}}
            list9mers[[i]] <- table9mers
          }else{  
            
            #codon=5           
            if(i==5){
              table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
              colnames(table9mers) <- AAseq$V1
              for (a in 1:i){       
                for (b in 1:ncol(table9mers)){
                  if(a==i-4){
                    table9mers[a,b] <- paste(colnames(table9mers)[b]
                                             , protseq[i+1]
                                             , protseq[i+2]
                                             , protseq[i+3]
                                             , protseq[i+4]
                                             , protseq[i+5]
                                             , protseq[i+6]
                                             , protseq[i+7]
                                             , protseq[i+8]
                                             ,sep="")
                  }else{if(a==i-3){
                    table9mers[a,b] <- paste(protseq[i-1]
                                             , colnames(table9mers)[b]
                                             , protseq[i+1]
                                             , protseq[i+2]
                                             , protseq[i+3]
                                             , protseq[i+4]
                                             , protseq[i+5]
                                             , protseq[i+6]
                                             , protseq[i+7]
                                             ,sep="")
                  }else{if(a==i-2){
                    table9mers[a,b] <- paste(protseq[i-2]
                                             , protseq[i-1]
                                             , colnames(table9mers)[b]
                                             , protseq[i+1]
                                             , protseq[i+2]
                                             , protseq[i+3]
                                             , protseq[i+4]
                                             , protseq[i+5]
                                             , protseq[i+6]
                                             ,sep="")
                  }else{if(a==i-1){
                    table9mers[a,b] <- paste(protseq[i-3]
                                             , protseq[i-2]
                                             , protseq[i-1]
                                             , colnames(table9mers)[b]
                                             , protseq[i+1]
                                             , protseq[i+2]
                                             , protseq[i+3]
                                             , protseq[i+4]
                                             , protseq[i+5]
                                             ,sep="")
                  }else{if(a==i){
                    table9mers[a,b] <- paste(protseq[i-4]
                                             , protseq[i-3]
                                             , protseq[i-2]
                                             , protseq[i-1]
                                             , colnames(table9mers)[b]
                                             , protseq[i+1]
                                             , protseq[i+2]
                                             , protseq[i+3]
                                             , protseq[i+4]
                                             ,sep="")}}}}}}}
              list9mers[[i]] <- table9mers
            }else{
              
              #codon=6           
              if(i==6){
                table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
                colnames(table9mers) <- AAseq$V1
                for (a in 1:i){       
                  for (b in 1:ncol(table9mers)){
                    if(a==i-5){
                      table9mers[a,b] <- paste(colnames(table9mers)[b]
                                               , protseq[i+1]
                                               , protseq[i+2]
                                               , protseq[i+3]
                                               , protseq[i+4]
                                               , protseq[i+5]
                                               , protseq[i+6]
                                               , protseq[i+7]
                                               , protseq[i+8]
                                               ,sep="")
                    }else{if(a==i-4){
                      table9mers[a,b] <- paste(protseq[i-1]
                                               , colnames(table9mers)[b]
                                               , protseq[i+1]
                                               , protseq[i+2]
                                               , protseq[i+3]
                                               , protseq[i+4]
                                               , protseq[i+5]
                                               , protseq[i+6]
                                               , protseq[i+7]
                                               ,sep="")
                    }else{if(a==i-3){
                      table9mers[a,b] <- paste(protseq[i-2]
                                               , protseq[i-1]
                                               , colnames(table9mers)[b]
                                               , protseq[i+1]
                                               , protseq[i+2]
                                               , protseq[i+3]
                                               , protseq[i+4]
                                               , protseq[i+5]
                                               , protseq[i+6]
                                               ,sep="")
                    }else{if(a==i-2){
                      table9mers[a,b] <- paste(protseq[i-3]
                                               , protseq[i-2]
                                               , protseq[i-1]
                                               , colnames(table9mers)[b]
                                               , protseq[i+1]
                                               , protseq[i+2]
                                               , protseq[i+3]
                                               , protseq[i+4]
                                               , protseq[i+5]
                                               ,sep="")
                    }else{if(a==i-1){
                      table9mers[a,b] <- paste(protseq[i-4]
                                               , protseq[i-3]
                                               , protseq[i-2]
                                               , protseq[i-1]
                                               , colnames(table9mers)[b]
                                               , protseq[i+1]
                                               , protseq[i+2]
                                               , protseq[i+3]
                                               , protseq[i+4]
                                               ,sep="")
                    }else{if(a==i){
                      table9mers[a,b] <- paste(protseq[i-5]
                                               , protseq[i-4]
                                               , protseq[i-3]
                                               , protseq[i-2]
                                               , protseq[i-1]
                                               , colnames(table9mers)[b]
                                               , protseq[i+1]
                                               , protseq[i+2]
                                               , protseq[i+3]
                                               ,sep="")}}}}}}}}
                list9mers[[i]] <- table9mers
              }else{    
                
                #codon=7           
                if(i==7){
                  table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
                  colnames(table9mers) <- AAseq$V1
                  for (a in 1:i){       
                    for (b in 1:ncol(table9mers)){
                      if(a==i-6){
                        table9mers[a,b] <- paste(colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 , protseq[i+3]
                                                 , protseq[i+4]
                                                 , protseq[i+5]
                                                 , protseq[i+6]
                                                 , protseq[i+7]
                                                 , protseq[i+8]
                                                 ,sep="")
                      }else{if(a==i-5){
                        table9mers[a,b] <- paste(protseq[i-1]
                                                 , colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 , protseq[i+3]
                                                 , protseq[i+4]
                                                 , protseq[i+5]
                                                 , protseq[i+6]
                                                 , protseq[i+7]
                                                 ,sep="")
                      }else{if(a==i-4){
                        table9mers[a,b] <- paste(protseq[i-2]
                                                 , protseq[i-1]
                                                 , colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 , protseq[i+3]
                                                 , protseq[i+4]
                                                 , protseq[i+5]
                                                 , protseq[i+6]
                                                 ,sep="")
                      }else{if(a==i-3){
                        table9mers[a,b] <- paste(protseq[i-3]
                                                 , protseq[i-2]
                                                 , protseq[i-1]
                                                 , colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 , protseq[i+3]
                                                 , protseq[i+4]
                                                 , protseq[i+5]
                                                 ,sep="")
                      }else{if(a==i-2){
                        table9mers[a,b] <- paste(protseq[i-4]
                                                 , protseq[i-3]
                                                 , protseq[i-2]
                                                 , protseq[i-1]
                                                 , colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 , protseq[i+3]
                                                 , protseq[i+4]
                                                 ,sep="")
                      }else{if(a==i-1){
                        table9mers[a,b] <- paste(protseq[i-5]
                                                 , protseq[i-4]
                                                 , protseq[i-3]
                                                 , protseq[i-2]
                                                 , protseq[i-1]
                                                 , colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 , protseq[i+3]
                                                 ,sep="")
                      }else{if(a==i){
                        table9mers[a,b] <- paste(protseq[i-6]
                                                 , protseq[i-5]
                                                 , protseq[i-4]
                                                 , protseq[i-3]
                                                 , protseq[i-2]
                                                 , protseq[i-1]
                                                 , colnames(table9mers)[b]
                                                 , protseq[i+1]
                                                 , protseq[i+2]
                                                 ,sep="")}}}}}}}}}
                  list9mers[[i]] <- table9mers
                }else{    
                  
                  #codon=8
                  if(i==8){
                    table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
                    colnames(table9mers) <- AAseq$V1
                    for (a in 1:i){       
                      for (b in 1:ncol(table9mers)){
                        if(a==i-7){
                          table9mers[a,b] <- paste(colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   , protseq[i+3]
                                                   , protseq[i+4]
                                                   , protseq[i+5]
                                                   , protseq[i+6]
                                                   , protseq[i+7]
                                                   , protseq[i+8]
                                                   ,sep="")
                        }else{if(a==i-6){
                          table9mers[a,b] <- paste(protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   , protseq[i+3]
                                                   , protseq[i+4]
                                                   , protseq[i+5]
                                                   , protseq[i+6]
                                                   , protseq[i+7]
                                                   ,sep="")
                        }else{if(a==i-5){
                          table9mers[a,b] <- paste(protseq[i-2]
                                                   , protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   , protseq[i+3]
                                                   , protseq[i+4]
                                                   , protseq[i+5]
                                                   , protseq[i+6]
                                                   ,sep="")
                        }else{if(a==i-4){
                          table9mers[a,b] <- paste(protseq[i-3]
                                                   , protseq[i-2]
                                                   , protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   , protseq[i+3]
                                                   , protseq[i+4]
                                                   , protseq[i+5]
                                                   ,sep="")
                        }else{if(a==i-3){
                          table9mers[a,b] <- paste(protseq[i-4]
                                                   , protseq[i-3]
                                                   , protseq[i-2]
                                                   , protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   , protseq[i+3]
                                                   , protseq[i+4]
                                                   ,sep="")
                        }else{if(a==i-2){
                          table9mers[a,b] <- paste(protseq[i-5]
                                                   , protseq[i-4]
                                                   , protseq[i-3]
                                                   , protseq[i-2]
                                                   , protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   , protseq[i+3]
                                                   ,sep="")
                        }else{if(a==i-1){
                          table9mers[a,b] <- paste(protseq[i-6]
                                                   , protseq[i-5]
                                                   , protseq[i-4]
                                                   , protseq[i-3]
                                                   , protseq[i-2]
                                                   , protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   , protseq[i+2]
                                                   ,sep="")
                        }else{if(a==i){
                          table9mers[a,b] <- paste(protseq[i-7]
                                                   , protseq[i-6]
                                                   , protseq[i-5]
                                                   , protseq[i-4]
                                                   , protseq[i-3]
                                                   , protseq[i-2]
                                                   , protseq[i-1]
                                                   , colnames(table9mers)[b]
                                                   , protseq[i+1]
                                                   ,sep="")}}}}}}}}}}
                    list9mers[[i]] <- table9mers
                  }else{   
                    
                    #codon=9
                    if(i==9){
                      table9mers <- as.data.frame(matrix(,nrow=i, ncol=AAn))
                      colnames(table9mers) <- AAseq$V1
                      for (a in 1:i){       
                        for (b in 1:ncol(table9mers)){
                          if(a==i-8){
                            table9mers[a,b] <- paste(colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     , protseq[i+3]
                                                     , protseq[i+4]
                                                     , protseq[i+5]
                                                     , protseq[i+6]
                                                     , protseq[i+7]
                                                     , protseq[i+8]
                                                     ,sep="")
                          }else{if(a==i-7){
                            table9mers[a,b] <- paste(protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     , protseq[i+3]
                                                     , protseq[i+4]
                                                     , protseq[i+5]
                                                     , protseq[i+6]
                                                     , protseq[i+7]
                                                     ,sep="")
                          }else{if(a==i-6){
                            table9mers[a,b] <- paste(protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     , protseq[i+3]
                                                     , protseq[i+4]
                                                     , protseq[i+5]
                                                     , protseq[i+6]
                                                     ,sep="")
                          }else{if(a==i-5){
                            table9mers[a,b] <- paste(protseq[i-3]
                                                     , protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     , protseq[i+3]
                                                     , protseq[i+4]
                                                     , protseq[i+5]
                                                     ,sep="")
                          }else{if(a==i-4){
                            table9mers[a,b] <- paste(protseq[i-4]
                                                     , protseq[i-3]
                                                     , protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     , protseq[i+3]
                                                     , protseq[i+4]
                                                     ,sep="")
                          }else{if(a==i-3){
                            table9mers[a,b] <- paste(protseq[i-5]
                                                     , protseq[i-4]
                                                     , protseq[i-3]
                                                     , protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     , protseq[i+3]
                                                     ,sep="")
                          }else{if(a==i-2){
                            table9mers[a,b] <- paste(protseq[i-6]
                                                     , protseq[i-5]
                                                     , protseq[i-4]
                                                     , protseq[i-3]
                                                     , protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     , protseq[i+2]
                                                     ,sep="")
                          }else{if(a==i-1){
                            table9mers[a,b] <- paste(protseq[i-7]
                                                     , protseq[i-6]
                                                     , protseq[i-5]
                                                     , protseq[i-4]
                                                     , protseq[i-3]
                                                     , protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     , protseq[i+1]
                                                     ,sep="")
                          }else{if(a==i){
                            table9mers[a,b] <- paste(protseq[i-8]
                                                     , protseq[i-7]
                                                     , protseq[i-6]
                                                     , protseq[i-5]
                                                     , protseq[i-4]
                                                     , protseq[i-3]
                                                     , protseq[i-2]
                                                     , protseq[i-1]
                                                     , colnames(table9mers)[b]
                                                     ,sep="")}}}}}}}}}}}
                      list9mers[[i]] <- table9mers
                    }else{       
                      
                      #codon=10 to codon=codon-9 
                      if(is.na(match(i,c(10:(codon-9)))>0)==F){
                        table9mers <- as.data.frame(matrix(,nrow=9, ncol=AAn))
                        colnames(table9mers) <- AAseq$V1
                        for (a in (i-8):i){       
                          for (b in 1:ncol(table9mers)){
                            if(a==i-8){
                              table9mers[1,b] <- paste(colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       , protseq[i+3]
                                                       , protseq[i+4]
                                                       , protseq[i+5]
                                                       , protseq[i+6]
                                                       , protseq[i+7]
                                                       , protseq[i+8]
                                                       ,sep="")
                            }else{if(a==i-7){
                              table9mers[2,b] <- paste(protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       , protseq[i+3]
                                                       , protseq[i+4]
                                                       , protseq[i+5]
                                                       , protseq[i+6]
                                                       , protseq[i+7]
                                                       ,sep="")
                            }else{if(a==i-6){
                              table9mers[3,b] <- paste(protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       , protseq[i+3]
                                                       , protseq[i+4]
                                                       , protseq[i+5]
                                                       , protseq[i+6]
                                                       ,sep="")
                            }else{if(a==i-5){
                              table9mers[4,b] <- paste(protseq[i-3]
                                                       , protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       , protseq[i+3]
                                                       , protseq[i+4]
                                                       , protseq[i+5]
                                                       ,sep="")
                            }else{if(a==i-4){
                              table9mers[5,b] <- paste(protseq[i-4]
                                                       , protseq[i-3]
                                                       , protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       , protseq[i+3]
                                                       , protseq[i+4]
                                                       ,sep="")
                            }else{if(a==i-3){
                              table9mers[6,b] <- paste(protseq[i-5]
                                                       , protseq[i-4]
                                                       , protseq[i-3]
                                                       , protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       , protseq[i+3]
                                                       ,sep="")
                            }else{if(a==i-2){
                              table9mers[7,b] <- paste(protseq[i-6]
                                                       , protseq[i-5]
                                                       , protseq[i-4]
                                                       , protseq[i-3]
                                                       , protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       , protseq[i+2]
                                                       ,sep="")
                            }else{if(a==i-1){
                              table9mers[8,b] <- paste(protseq[i-7]
                                                       , protseq[i-6]
                                                       , protseq[i-5]
                                                       , protseq[i-4]
                                                       , protseq[i-3]
                                                       , protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       , protseq[i+1]
                                                       ,sep="")
                            }else{if(a==i){
                              table9mers[9,b] <- paste(protseq[i-8]
                                                       , protseq[i-7]
                                                       , protseq[i-6]
                                                       , protseq[i-5]
                                                       , protseq[i-4]
                                                       , protseq[i-3]
                                                       , protseq[i-2]
                                                       , protseq[i-1]
                                                       , colnames(table9mers)[b]
                                                       ,sep="")}}}}}}}}}}}
                        list9mers[[i]] <- table9mers
                      }else{       
                        
                        #codon=codon-8
                        if(i==codon-8){
                          table9mers <- as.data.frame(matrix(,nrow=9, ncol=AAn))
                          colnames(table9mers) <- AAseq$V1
                          for (a in i:codon){       
                            for (b in 1:ncol(table9mers)){
                              if(a==i+8){
                                table9mers[9,b] <- paste(protseq[i-8]
                                                         , protseq[i-7]
                                                         , protseq[i-6]
                                                         , protseq[i-5]
                                                         , protseq[i-4]
                                                         , protseq[i-3]
                                                         , protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         ,sep="")
                              }else{if(a==i+7){
                                table9mers[8,b] <- paste(protseq[i-7]
                                                         , protseq[i-6]
                                                         , protseq[i-5]
                                                         , protseq[i-4]
                                                         , protseq[i-3]
                                                         , protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         ,sep="")
                              }else{if(a==i+6){
                                table9mers[7,b] <- paste(protseq[i-6]
                                                         , protseq[i-5]
                                                         , protseq[i-4]
                                                         , protseq[i-3]
                                                         , protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         ,sep="")
                              }else{if(a==i+5){
                                table9mers[6,b] <- paste(protseq[i-5]
                                                         , protseq[i-4]
                                                         , protseq[i-3]
                                                         , protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         , protseq[i+3]
                                                         ,sep="")
                              }else{if(a==i+4){
                                table9mers[5,b] <- paste(protseq[i-4]
                                                         , protseq[i-3]
                                                         , protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         , protseq[i+3]
                                                         , protseq[i+4]
                                                         ,sep="")
                              }else{if(a==i+3){
                                table9mers[4,b] <- paste(protseq[i-3]
                                                         , protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         , protseq[i+3]
                                                         , protseq[i+4]
                                                         , protseq[i+5]
                                                         ,sep="")
                              }else{if(a==i+2){
                                table9mers[3,b] <- paste(protseq[i-2]
                                                         , protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         , protseq[i+3]
                                                         , protseq[i+4]
                                                         , protseq[i+5]
                                                         , protseq[i+6]
                                                         ,sep="")
                              }else{if(a==i+1){
                                table9mers[2,b] <- paste(protseq[i-1]
                                                         , colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         , protseq[i+3]
                                                         , protseq[i+4]
                                                         , protseq[i+5]
                                                         , protseq[i+6]
                                                         , protseq[i+7]
                                                         ,sep="")
                              }else{if(a==i){
                                table9mers[1,b] <- paste(colnames(table9mers)[b]
                                                         , protseq[i+1]
                                                         , protseq[i+2]
                                                         , protseq[i+3]
                                                         , protseq[i+4]
                                                         , protseq[i+5]
                                                         , protseq[i+6]
                                                         , protseq[i+7]
                                                         , protseq[i+8]
                                                         ,sep="")}}}}}}}}}}}
                          list9mers[[i]] <- table9mers   
                        }else{  
                          
                          #codon=codon-7
                          if(i==codon-7){
                            table9mers <- as.data.frame(matrix(,nrow=8, ncol=AAn))
                            colnames(table9mers) <- AAseq$V1
                            for (a in i:codon){       
                              for (b in 1:ncol(table9mers)){
                                if(a==i+7){
                                  table9mers[8,b] <- paste(protseq[i-8]
                                                           , protseq[i-7]
                                                           , protseq[i-6]
                                                           , protseq[i-5]
                                                           , protseq[i-4]
                                                           , protseq[i-3]
                                                           , protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           ,sep="")
                                }else{if(a==i+6){
                                  table9mers[7,b] <- paste(protseq[i-7]
                                                           , protseq[i-6]
                                                           , protseq[i-5]
                                                           , protseq[i-4]
                                                           , protseq[i-3]
                                                           , protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           ,sep="")
                                }else{if(a==i+5){
                                  table9mers[6,b] <- paste(protseq[i-6]
                                                           , protseq[i-5]
                                                           , protseq[i-4]
                                                           , protseq[i-3]
                                                           , protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           , protseq[i+2]
                                                           ,sep="")
                                }else{if(a==i+4){
                                  table9mers[5,b] <- paste(protseq[i-5]
                                                           , protseq[i-4]
                                                           , protseq[i-3]
                                                           , protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           , protseq[i+2]
                                                           , protseq[i+3]
                                                           ,sep="")
                                }else{if(a==i+3){
                                  table9mers[4,b] <- paste(protseq[i-4]
                                                           , protseq[i-3]
                                                           , protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           , protseq[i+2]
                                                           , protseq[i+3]
                                                           , protseq[i+4]
                                                           ,sep="")
                                }else{if(a==i+2){
                                  table9mers[3,b] <- paste(protseq[i-3]
                                                           , protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           , protseq[i+2]
                                                           , protseq[i+3]
                                                           , protseq[i+4]
                                                           , protseq[i+5]
                                                           ,sep="")
                                }else{if(a==i+1){
                                  table9mers[2,b] <- paste(protseq[i-2]
                                                           , protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           , protseq[i+2]
                                                           , protseq[i+3]
                                                           , protseq[i+4]
                                                           , protseq[i+5]
                                                           , protseq[i+6]
                                                           ,sep="")
                                }else{if(a==i){
                                  table9mers[1,b] <- paste(protseq[i-1]
                                                           , colnames(table9mers)[b]
                                                           , protseq[i+1]
                                                           , protseq[i+2]
                                                           , protseq[i+3]
                                                           , protseq[i+4]
                                                           , protseq[i+5]
                                                           , protseq[i+6]
                                                           , protseq[i+7]
                                                           ,sep="")}}}}}}}}}}
                            list9mers[[i]] <- table9mers   
                          }else{  
                            
                            #codon=codon-6
                            if(i==codon-6){
                              table9mers <- as.data.frame(matrix(,nrow=7, ncol=AAn))
                              colnames(table9mers) <- AAseq$V1
                              for (a in i:codon){       
                                for (b in 1:ncol(table9mers)){
                                  if(a==i+6){
                                    table9mers[7,b] <- paste(protseq[i-8]
                                                             , protseq[i-7]
                                                             , protseq[i-6]
                                                             , protseq[i-5]
                                                             , protseq[i-4]
                                                             , protseq[i-3]
                                                             , protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             ,sep="")
                                  }else{if(a==i+5){
                                    table9mers[6,b] <- paste(protseq[i-7]
                                                             , protseq[i-6]
                                                             , protseq[i-5]
                                                             , protseq[i-4]
                                                             , protseq[i-3]
                                                             , protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             , protseq[i+1]
                                                             ,sep="")
                                  }else{if(a==i+4){
                                    table9mers[5,b] <- paste(protseq[i-6]
                                                             , protseq[i-5]
                                                             , protseq[i-4]
                                                             , protseq[i-3]
                                                             , protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             , protseq[i+1]
                                                             , protseq[i+2]
                                                             ,sep="")
                                  }else{if(a==i+3){
                                    table9mers[4,b] <- paste(protseq[i-5]
                                                             , protseq[i-4]
                                                             , protseq[i-3]
                                                             , protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             , protseq[i+1]
                                                             , protseq[i+2]
                                                             , protseq[i+3]
                                                             ,sep="")
                                  }else{if(a==i+2){
                                    table9mers[3,b] <- paste(protseq[i-4]
                                                             , protseq[i-3]
                                                             , protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             , protseq[i+1]
                                                             , protseq[i+2]
                                                             , protseq[i+3]
                                                             , protseq[i+4]
                                                             ,sep="")
                                  }else{if(a==i+1){
                                    table9mers[2,b] <- paste(protseq[i-3]
                                                             , protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             , protseq[i+1]
                                                             , protseq[i+2]
                                                             , protseq[i+3]
                                                             , protseq[i+4]
                                                             , protseq[i+5]
                                                             ,sep="")
                                  }else{if(a==i){
                                    table9mers[1,b] <- paste(protseq[i-2]
                                                             , protseq[i-1]
                                                             , colnames(table9mers)[b]
                                                             , protseq[i+1]
                                                             , protseq[i+2]
                                                             , protseq[i+3]
                                                             , protseq[i+4]
                                                             , protseq[i+5]
                                                             , protseq[i+6]
                                                             ,sep="")}}}}}}}}}
                              list9mers[[i]] <- table9mers   
                            }else{  
                              
                              #codon=codon-5
                              if(i==codon-5){
                                table9mers <- as.data.frame(matrix(,nrow=6, ncol=AAn))
                                colnames(table9mers) <- AAseq$V1
                                for (a in i:codon){       
                                  for (b in 1:ncol(table9mers)){
                                    if(a==i+5){
                                      table9mers[6,b] <- paste(protseq[i-8]
                                                               , protseq[i-7]
                                                               , protseq[i-6]
                                                               , protseq[i-5]
                                                               , protseq[i-4]
                                                               , protseq[i-3]
                                                               , protseq[i-2]
                                                               , protseq[i-1]
                                                               , colnames(table9mers)[b]
                                                               ,sep="")
                                    }else{if(a==i+4){
                                      table9mers[5,b] <- paste(protseq[i-7]
                                                               , protseq[i-6]
                                                               , protseq[i-5]
                                                               , protseq[i-4]
                                                               , protseq[i-3]
                                                               , protseq[i-2]
                                                               , protseq[i-1]
                                                               , colnames(table9mers)[b]
                                                               , protseq[i+1]
                                                               ,sep="")
                                    }else{if(a==i+3){
                                      table9mers[4,b] <- paste(protseq[i-6]
                                                               , protseq[i-5]
                                                               , protseq[i-4]
                                                               , protseq[i-3]
                                                               , protseq[i-2]
                                                               , protseq[i-1]
                                                               , colnames(table9mers)[b]
                                                               , protseq[i+1]
                                                               , protseq[i+2]
                                                               ,sep="")
                                    }else{if(a==i+2){
                                      table9mers[3,b] <- paste(protseq[i-5]
                                                               , protseq[i-4]
                                                               , protseq[i-3]
                                                               , protseq[i-2]
                                                               , protseq[i-1]
                                                               , colnames(table9mers)[b]
                                                               , protseq[i+1]
                                                               , protseq[i+2]
                                                               , protseq[i+3]
                                                               ,sep="")
                                    }else{if(a==i+1){
                                      table9mers[2,b] <- paste(protseq[i-4]
                                                               , protseq[i-3]
                                                               , protseq[i-2]
                                                               , protseq[i-1]
                                                               , colnames(table9mers)[b]
                                                               , protseq[i+1]
                                                               , protseq[i+2]
                                                               , protseq[i+3]
                                                               , protseq[i+4]
                                                               ,sep="")
                                    }else{if(a==i){
                                      table9mers[1,b] <- paste(protseq[i-3]
                                                               , protseq[i-2]
                                                               , protseq[i-1]
                                                               , colnames(table9mers)[b]
                                                               , protseq[i+1]
                                                               , protseq[i+2]
                                                               , protseq[i+3]
                                                               , protseq[i+4]
                                                               , protseq[i+5]
                                                               ,sep="")}}}}}}}}
                                list9mers[[i]] <- table9mers   
                              }else{ 
                                
                                #codon=codon-4
                                if(i==codon-4){
                                  table9mers <- as.data.frame(matrix(,nrow=5, ncol=AAn))
                                  colnames(table9mers) <- AAseq$V1
                                  for (a in i:codon){       
                                    for (b in 1:ncol(table9mers)){
                                      if(a==i+4){
                                        table9mers[5,b] <- paste(protseq[i-8]
                                                                 , protseq[i-7]
                                                                 , protseq[i-6]
                                                                 , protseq[i-5]
                                                                 , protseq[i-4]
                                                                 , protseq[i-3]
                                                                 , protseq[i-2]
                                                                 , protseq[i-1]
                                                                 , colnames(table9mers)[b]
                                                                 ,sep="")
                                      }else{if(a==i+3){
                                        table9mers[4,b] <- paste(protseq[i-7]
                                                                 , protseq[i-6]
                                                                 , protseq[i-5]
                                                                 , protseq[i-4]
                                                                 , protseq[i-3]
                                                                 , protseq[i-2]
                                                                 , protseq[i-1]
                                                                 , colnames(table9mers)[b]
                                                                 , protseq[i+1]
                                                                 ,sep="")
                                      }else{if(a==i+2){
                                        table9mers[3,b] <- paste(protseq[i-6]
                                                                 , protseq[i-5]
                                                                 , protseq[i-4]
                                                                 , protseq[i-3]
                                                                 , protseq[i-2]
                                                                 , protseq[i-1]
                                                                 , colnames(table9mers)[b]
                                                                 , protseq[i+1]
                                                                 , protseq[i+2]
                                                                 ,sep="")
                                      }else{if(a==i+1){
                                        table9mers[2,b] <- paste(protseq[i-5]
                                                                 , protseq[i-4]
                                                                 , protseq[i-3]
                                                                 , protseq[i-2]
                                                                 , protseq[i-1]
                                                                 , colnames(table9mers)[b]
                                                                 , protseq[i+1]
                                                                 , protseq[i+2]
                                                                 , protseq[i+3]
                                                                 ,sep="")
                                      }else{if(a==i){
                                        table9mers[1,b] <- paste(protseq[i-4]
                                                                 , protseq[i-3]
                                                                 , protseq[i-2]
                                                                 , protseq[i-1]
                                                                 , colnames(table9mers)[b]
                                                                 , protseq[i+1]
                                                                 , protseq[i+2]
                                                                 , protseq[i+3]
                                                                 , protseq[i+4]
                                                                 ,sep="")}}}}}}}
                                  list9mers[[i]] <- table9mers   
                                }else{
                                  
                                  #codon=codon-3
                                  if(i==codon-3){
                                    table9mers <- as.data.frame(matrix(,nrow=4, ncol=AAn))
                                    colnames(table9mers) <- AAseq$V1
                                    for (a in i:codon){       
                                      for (b in 1:ncol(table9mers)){
                                        if(a==i+3){
                                          table9mers[4,b] <- paste(protseq[i-8]
                                                                   , protseq[i-7]
                                                                   , protseq[i-6]
                                                                   , protseq[i-5]
                                                                   , protseq[i-4]
                                                                   , protseq[i-3]
                                                                   , protseq[i-2]
                                                                   , protseq[i-1]
                                                                   , colnames(table9mers)[b]
                                                                   ,sep="")
                                        }else{if(a==i+2){
                                          table9mers[3,b] <- paste(protseq[i-7]
                                                                   , protseq[i-6]
                                                                   , protseq[i-5]
                                                                   , protseq[i-4]
                                                                   , protseq[i-3]
                                                                   , protseq[i-2]
                                                                   , protseq[i-1]
                                                                   , colnames(table9mers)[b]
                                                                   , protseq[i+1]
                                                                   ,sep="")
                                        }else{if(a==i+1){
                                          table9mers[2,b] <- paste(protseq[i-6]
                                                                   , protseq[i-5]
                                                                   , protseq[i-4]
                                                                   , protseq[i-3]
                                                                   , protseq[i-2]
                                                                   , protseq[i-1]
                                                                   , colnames(table9mers)[b]
                                                                   , protseq[i+1]
                                                                   , protseq[i+2]
                                                                   ,sep="")
                                        }else{if(a==i){
                                          table9mers[1,b] <- paste(protseq[i-5]
                                                                   , protseq[i-4]
                                                                   , protseq[i-3]
                                                                   , protseq[i-2]
                                                                   , protseq[i-1]
                                                                   , colnames(table9mers)[b]
                                                                   , protseq[i+1]
                                                                   , protseq[i+2]
                                                                   , protseq[i+3]
                                                                   ,sep="")}}}}}}
                                    list9mers[[i]] <- table9mers   
                                  }else{
                                    
                                    #codon=codon-2
                                    if(i==codon-2){
                                      table9mers <- as.data.frame(matrix(,nrow=3, ncol=AAn))
                                      colnames(table9mers) <- AAseq$V1
                                      for (a in i:codon){       
                                        for (b in 1:ncol(table9mers)){
                                          if(a==i+2){
                                            table9mers[3,b] <- paste(protseq[i-8]
                                                                     , protseq[i-7]
                                                                     , protseq[i-6]
                                                                     , protseq[i-5]
                                                                     , protseq[i-4]
                                                                     , protseq[i-3]
                                                                     , protseq[i-2]
                                                                     , protseq[i-1]
                                                                     , colnames(table9mers)[b]
                                                                     ,sep="")
                                          }else{if(a==i+1){
                                            table9mers[2,b] <- paste(protseq[i-7]
                                                                     , protseq[i-6]
                                                                     , protseq[i-5]
                                                                     , protseq[i-4]
                                                                     , protseq[i-3]
                                                                     , protseq[i-2]
                                                                     , protseq[i-1]
                                                                     , colnames(table9mers)[b]
                                                                     , protseq[i+1]
                                                                     ,sep="")
                                          }else{if(a==i){
                                            table9mers[1,b] <- paste(protseq[i-6]
                                                                     , protseq[i-5]
                                                                     , protseq[i-4]
                                                                     , protseq[i-3]
                                                                     , protseq[i-2]
                                                                     , protseq[i-1]
                                                                     , colnames(table9mers)[b]
                                                                     , protseq[i+1]
                                                                     , protseq[i+2]
                                                                     ,sep="")}}}}}
                                      list9mers[[i]] <- table9mers   
                                    }else{
                                      
                                      #codon=codon-1
                                      if(i==codon-1){
                                        table9mers <- as.data.frame(matrix(,nrow=2, ncol=AAn))
                                        colnames(table9mers) <- AAseq$V1
                                        for (a in i:codon){       
                                          for (b in 1:ncol(table9mers)){
                                            if(a==i+1){
                                              table9mers[2,b] <- paste(protseq[i-8]
                                                                       , protseq[i-7]
                                                                       , protseq[i-6]
                                                                       , protseq[i-5]
                                                                       , protseq[i-4]
                                                                       , protseq[i-3]
                                                                       , protseq[i-2]
                                                                       , protseq[i-1]
                                                                       , colnames(table9mers)[b]
                                                                       ,sep="")
                                            }else{if(a==i){
                                              table9mers[1,b] <- paste(protseq[i-7]
                                                                       , protseq[i-6]
                                                                       , protseq[i-5]
                                                                       , protseq[i-4]
                                                                       , protseq[i-3]
                                                                       , protseq[i-2]
                                                                       , protseq[i-1]
                                                                       , colnames(table9mers)[b]
                                                                       , protseq[i+1]
                                                                       ,sep="")}}}}
                                        list9mers[[i]] <- table9mers   
                                      }else{
                                        #codon=codon
                                        if(i==codon){
                                          table9mers <- as.data.frame(matrix(,nrow=1, ncol=AAn))
                                          colnames(table9mers) <- AAseq$V1
                                          
                                          for (a in i){ 
                                            for (b in 1:ncol(table9mers)){
                                              table9mers[1,b] <- paste(protseq[i-8]
                                                                       , protseq[i-7]
                                                                       , protseq[i-6]
                                                                       , protseq[i-5]
                                                                       , protseq[i-4]
                                                                       , protseq[i-3]
                                                                       , protseq[i-2]
                                                                       , protseq[i-1]
                                                                       , colnames(table9mers)[b]
                                                                       ,sep="")}}
                                          list9mers[[i]] <- table9mers
                                        }else{ }}}}}}}}}}}}}}}}}}}
    
    
    
    
    write.csv2(list9mers[[i]], file=paste("files/9mers_codon",i,".csv",sep="")) 

  }
}

  #----