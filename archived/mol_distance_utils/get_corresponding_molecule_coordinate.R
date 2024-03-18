getqrycoor<-function(rcmap, xmap, labelstart){

# source("/home/users/elam/rscripts/readmaps.R")
# source("/home/users/elam/rscripts/paste0.R")
# r<-readcmap(rcmap)
# x<-readxmap(xmap)
  r<-rcmap
  x<-xmap

refstart<-r[r$SiteID == labelstart,]
# refstartpos<-refstart$Position
refstartpos<-data.frame(as.numeric(as.character(refstart$CMapId)), refstart$Position)
refstartmolcount<-refstart$Occurrence

# refend<-r[r$SiteID == labelend,]
# refendpos<-refend$Position
# refendmolcount<-refend$Occurrence

withalignments<-x[grepl(labelstart, x$Alignment),]
#withalignmentsend<-x[grepl(labelend, x$Alignment),]

correspondingmolcoor<-function(rownum, reflabelnum, withalignments_1){
  print(rownum)
  permolalignment<-withalignments_1[rownum,]
  eachlabel<-strsplit(permolalignment$Alignment, split="[)(]")
  wantedqrylabels<-eachlabel
  #wantedqrylabels<-eachlabel[grepl(reflabelnum, eachlabel)]
  #print(length(wantedqrylabels))
  
  getperwantedqrylabel<-function(num, reflabelnum){
    #print(num)
    perwantedqrylabels<-wantedqrylabels[[num]]
    #print("pass1")
    the_align<-perwantedqrylabels[grep(reflabelnum, perwantedqrylabels)]
    the_align<-unlist(strsplit(the_align, split=". "), recursive=F)
    peralign<-strsplit(the_align, split=",")
    findcorrespond<-data.frame(matrix(unlist(peralign), length(peralign), byrow =T))
    good<-which(findcorrespond$X1 == reflabelnum)
    #peralign<-data.frame(matrix(unlist(peralign), nrow=length(peralign), byrow=T))
    if(length(good)>0){
      the_align<-the_align[good]
      qrylabel<-do.call(rbind, strsplit(the_align, split=","))[,2]
      #wantedqrylabel<-data.frame(strsplit(perwantedqrylabels, ","))
    
      print(paste0(length(qrylabel), " found labels"))
      #print("pass2")
      return(qrylabel)
    } else {
      return("NA")
    }
  }
  #print(length(wantedqrylabels))
  qrylabels<-lapply(seq(1, length(wantedqrylabels), 1), getperwantedqrylabel, reflabelnum)
  # print("qrylabels done")
  return(qrylabels)
}
  num_row<-as.integer(nrow(withalignments))  
  qrycoor<-lapply(seq(1, num_row, 1), correspondingmolcoor, labelstart, withalignments_1=withalignments)
  print("qrycoor")

  # qrycoor<-lapply(seq(1, nrow(withalignments), 1), correspondingmolcoor, labelstart)
  # qryends<-lapply(seq(1, nrow(withalignmentsend), 1), correspondingmolcoor, refend, withalignmentsend)
  #print(paste0(nrow(withalignments), " alignments"))
  #print(paste0(length(qrycoor)," qry coordinates"))


  # result<-data.frame(withalignments$QryContigID, do.call(rbind, qrycoor))
  result<-data.frame(withalignments$RefContigID ,withalignments$QryContigID, do.call(rbind, qrycoor))
  colnames(result)<-c("RefContigID", "QryContigID", "Qrylabel")
  print("result data frame initiated")
  # result$refPos<-refstartpos
  # result$Qrylabel<-unlist(result$Qrylabel)
  result$refPos<-refstartpos$Position[match(result$RefContigID, refstartpos$CMapId)]
  result$Qrylabel<-unlist(result$Qrylabel)
  result<-result[!result$Qrylabel=="NA",]
  
  return(result)
}


