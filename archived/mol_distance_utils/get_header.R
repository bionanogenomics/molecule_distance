getcolnames<-function(rarepath){
  rare1<-read.table(rarepath, comment.char="#", sep="\t", stringsAsFactors = F)
  nhl<-system(paste0("grep '^#' ", rarepath, " | wc -l"), intern =TRUE)
  #header<-readLines(rarepath, n=nhl)
  header<-system(paste0("less ", rarepath, "| grep '^#h'"), intern=TRUE)
  if(length(header)==0){
    header<-system(paste0("less ", rarepath, "| grep '^# SVIndex'"), intern=TRUE)
  } 
  if(length(header)==0){
    header<-system(paste0("less ", rarepath, "| grep '^#clusterId'"), intern=TRUE)
  } 
  if(length(header)==0){
    header<-system(paste0("less ", rarepath, "| grep '^#'"), intern=TRUE)
  } 
  colnames(rare1)<-strsplit(header, split="\t")[[1]]
  #colnames(rare1)<-strsplit(header[as.integer(nhl)], split="\t")[[1]]
  return(rare1)
}

getheader<-function(rarepath){
  rare1<-read.table(rarepath, comment.char="#", sep="\t", stringsAsFactors = F)
  nhl<-system(paste0("grep '^#' ", rarepath, " | wc -l"), intern =TRUE)
  header<-readLines(rarepath, n=nhl)
  #colnames(rare1)<-strsplit(header[as.integer(nhl)], split="\t")[[1]]
  return(header)
}

outputsmap<-function(dataframe, headerfrom, suffix){
  header2<-getheader(headerfrom)
  prefix<-strsplit(headerfrom, split="[.]")[[1]][1]
  type<-strsplit(headerfrom, split="[.]")[[1]][2]
  cat(header2, file =paste0(prefix, "_", suffix,".",type), sep="\n")
  write.table(dataframe, file =paste0(prefix, "_", suffix,".",type), append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  
}


readtable<-function(filepath){
  see<-readLines(filepath)
  nhl<-system(paste0("grep '^#' ", filepath, " | wc -l"), intern =TRUE)
  rows<-seq(as.integer(nhl)+1, length(see),1)
  
  splitlines<-function(rownumber){
    perrow<-see[rownumber]
    eachrow<-unlist(strsplit(perrow, split = "\t"))
    return (eachrow)
  }
  
  table<-data.frame(do.call(rbind, lapply( rows,splitlines)))
  
  header<-system(paste0("less ", filepath, "| grep '^#h'"), intern=TRUE)
  if(length(header)==0){
    header<-system(paste0("less ", filepath, "| grep '^# SVIndex'"), intern=TRUE)
  } 
    
  colnames(table)<-strsplit(header, split="\t")[[1]]
  
  type<-system(paste0("less ", filepath, "| grep '^#f'"), intern=TRUE)
  types<-unlist(strsplit(type, split= "\t"))
  strings<-which(types == "string")
  stringheaders<-unlist(strsplit(header, split="\t"))[strings]
  int<-which(types != "string")
  intheader<-unlist(strsplit(header, split="\t"))[int]
  
  # table[, colnames(table) %in% stringheaders]<-as.character(table[,colnames(table) %in% stringheaders])
  # table[, colnames(table) %in% intheader]<-as.integer(as.character((table[,colnames(table) %in% intheader])))
  # 
  # 
  # bestrings<-function(n){
  #   colnum<-stringheaders[n]
  #   table[, colnames(table) == colnum]<-as.character(table[,colnames(table)==colnum])
  # }
  # 
  # beinteger<-function(intnum){
  #   table[,intnum]<-as.integer(as.character(table[,intnum]))
  # }
  # 
  # lapply(seq(1, length(stringheaders),1), bestrings)
  # lapply(int, beinteger)

return(table)
}
