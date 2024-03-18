#' Plot distance between labels of molecules that align to chosen reference label. 
#'
#' @name mol_distance_distribution
#' @keywords distance
#' @param rcmappath: anchor/reference cmap
#' @param xmap: xmap with alignment information
#' @param qcmappath: query/molecule cmap
#' @param startlabelid: ID of left reference label of interest
#' @param endlabelid: ID of right reference label of interest
#' @param binw: Binwidth for histogram. Default is 400
#'# @param kmax: Maximum number of clusters. Currently not executing this limit. 
#' @return A plot/histogram with distances in molecules between labels that align to the chosen reference labels. 
#' @export

#library(mixtools)
library(mclust)
library(prabclus)
library(ggplot2)
source("/home/users6/sshukor/molecule_distance_scripts/JLee_scripts/get_corresponding_molecule_coordinate.R")
source("/home/users6/sshukor/molecule_distance_scripts/JLee_scripts/get_header.R")

mol_distance_distribution_20210614_test_filtering<-function(rcmappath, xmap, qcmappath, startlabelid, endlabelid, nk=NULL, additional_noise_removal=FALSE, binw=400){
  
# print(paste0(rcmappath, "_", startlabelid))
qcmap<-getcolnames(qcmappath)
qcmap$CMapId<-as.character(qcmap$`#h CMapId`)

r<-getcolnames(rcmappath)
r$CMapId<-as.character(r$`#h CMapId`)
x<-read.table(xmap, stringsAsFactors=F)
h<-getheader(xmap)
h2<-h[length(h)-1]
h<-strsplit(h2, split="\t")[[1]]
colnames(x)<-h

starts<-getqrycoor(r, x, startlabelid)
ends<-getqrycoor(r, x, endlabelid)

###Downsample to the most 1000 molecules
if(dim(starts)[1]>1000){
  starts<-starts[sample(nrow(starts),1000),]
  ends<-ends[sample(nrow(ends), 1000),]
}

###temporary -get the righmost label if the mols that have multiple alignments to one ref label
starts$Qrylabel<-unlist(lapply(starts$Qrylabel, max))
###temporary - get the leftmost label if the mols have multiple alignments to one ref label
ends$Qrylabel<-unlist(lapply(ends$Qrylabel, min))

getcoor<-function(n, starts){
  print(n)
  perqry<-starts[n,]
  findcoor<-qcmap[qcmap$CMapId == perqry$QryContigID,]
  findcoor<-findcoor[findcoor$SiteID == perqry$Qrylabel,]$Position
  return(findcoor)
}

is.error <- function(x) inherits(x, "try-error")

###Get Qry position
if(nrow(starts)<10 | nrow(ends)<10){   ###Too few molecules at the start of end (Added on 04/20/2021)
  result<-c(r$CMapId[1], startlabelid, starts$refpos[1], min(nrow(starts), nrow(ends)), 0, 0, 0, 0)
  print(paste("Only <10 molecules are the label start or label end"))
  return (result)
}

starts$QryPos<-lapply(seq(1, nrow(starts), 1), getcoor, starts)
ends$QryPos<-lapply(seq(1, nrow(ends),1), getcoor, ends)
completeframe <- merge(starts,ends,by="QryContigID")
print(paste0(nrow(starts)," molecules at start, ", nrow(ends), " molecules at end, ", nrow (completeframe), " spanned both labels"))

###x is left breakpoint, y is right breakpoint
extramol<-which(completeframe$QryPos.x == "numeric(0)"|completeframe$QryPos.y == "numeric(0)")  #Remove any mol that is not in q.cmap (caused by potential RefAligner bug)
if(length(extramol)>0){
  completeframe<-completeframe[-extramol,]
}
completeframe$QryPos.x<-unlist(completeframe$QryPos.x)
completeframe$QryPos.y<-unlist(completeframe$QryPos.y)

if(nrow(completeframe)<10){   ###Too few molecules to be meaningful
  k<-0
  kmax<-0
  result<-c(r$CMapId[1], startlabelid, completeframe$refdistance[1], nrow(completeframe), k, rep(0, kmax), rep(0, kmax), rep(0, kmax))
  print(paste("Only", nrow(completeframe), " molecules span the labels"))
  return (result)
} else if(nrow(completeframe)>10){
  
  completeframe$distance<-abs(as.numeric(as.character(completeframe$QryPos.y))-as.numeric(as.character(completeframe$QryPos.x)))
  completeframe$refstartlabelid<-startlabelid
  completeframe$refendlabelid<-endlabelid
  completeframe$refdistance<-abs(as.numeric(as.character(completeframe$refpos.y))-as.numeric(as.character(completeframe$refpos.x)))
  
  completeframe$distance<-round(completeframe$distance,0)
  completeframeoriginal<-round(completeframe$distance,0)
  
  print("pass")
  
  ###Try to clean up data before clustering - Added 04/02/2021
  if(additional_noise_removal==TRUE){         ###Typically only needed if two clusters are close and have issue separating due to noise at further regions
    p1<-ggplot(data.frame(completeframeoriginal), aes(x=completeframeoriginal))+geom_histogram(aes(x=completeframeoriginal, y=..count..), fill=NA, color="black")
    pdc<-ggplot_build(p1)
    #pdc$data[[1]]
  
    cumulativecount<-cumsum(pdc$data[[1]]$count)
    removal=c()
    for (ci in seq(1,5,1)){
      if (cumulativecount[ci] == cumulativecount[ci+5] & pdc$data[[1]]$count[ci]==1){     ###checking left most 5 bins for non-zero bins separated by more than 5 zero bins (noise removal)
        removal<-c(removal, ci)
        completeframe<-completeframe[!(completeframe$distance >= pdc$data[[1]]$xmin[ci] & completeframe$distance <= pdc$data[[1]]$xmax[ci]),]
      }
    }
  
    for (ci in seq(length(cumulativecount),length(cumulativecount)-5,-1)){
      if (cumulativecount[ci-1] == cumulativecount[ci-5] & pdc$data[[1]]$count[ci]==1){     ###checking right most 5 bins for non-zero bins separated by more than 5 zero bins (noise removal)
        removal<-c(removal, ci)
        completeframe<-completeframe[!(completeframe$distance >= pdc$data[[1]]$xmin[ci] & completeframe$distance <= pdc$data[[1]]$xmax[ci]),]
      }
    }
    print(removal)
    print(paste(nrow(completeframe), " molecules left after noise filtering from checking the left most and right most 5 bins"))
  }
  ###
  
  mod <- Mclust(completeframe$distance, modelNames="V")
  if(!is.null(nk)){
    mod<-Mclust(completeframe$distance, modelNames="V", G=nk)
  }
  j<-summary(mod)
  summary(mod$BIC)
  mclust1Dplot(data = completeframe$distance, what = "density", parameters=mod$parameters, z=mod$z)
  abline(v=mod$parameters$mean, lty=3)
  mclust1Dplot(data = completeframe$distance, what = "classification", parameters=mod$parameters, z=mod$z)
  mclust1Dplot(data = completeframe$distance, what = "uncertainty", parameters=mod$parameters, z=mod$z)  ###uncertainty is how uncertain it is to be asscoiated with the group classified. 
  plot.mclustBIC(mod$BIC)
  k<-mod$G
  print(nrow(completeframe))
  print(paste(completeframe$refdistance[1], "ref distance"))
  print(mod$parameters$mean)
  print(mod$parameters$pro)
  print(mod$sd)
  print(mod$parameters$variance$sigmasq)
  hist(mod$uncertainty)
  print("done round 1")
  colors<-c("red", "green", "blue")
  
  
  ###Clear data in outlier cluster or reduce number of clusters
  dis<-data.frame(table(mod$classification))

###Plotting function
    Plotting_hist_and_Gau<-function(mod){
      newcolor<-array()
      print(colors)
      adjusteddnorm<-function(data, mm, sd, prop){
        re<-prop*dnorm(data, mean=mm, sd=sd)*nrow(completeframe)*binw
        return(re)
      }
      
      p1<-ggplot(data.frame(completeframeoriginal), aes(x=completeframeoriginal))+geom_histogram(aes(x=completeframeoriginal, y=..count..), fill=NA, color="black", binwidth=binw)
      
      print(length(mod$parameters$mean))
      
      ###Gather Component info for legend
      legwrite=array()
      for (k in 1:length(mod$parameters$mean)){
        print(k)
        newleg<-paste0("C-", nrow(completeframe), "\nD-", completeframe$refdistance[1], "\nM-", format(mod$parameters$mean[k],digits=3), "\nSD-", format(mod$sd[k], digits=3), "\nSD/M-", format((mod$sd[k]/mod$parameters$mean[k]*100), digit=3), "%\nP-", format(mod$parameters$pro[k], digits=2))
        print(newleg)
        legwrite<-append(legwrite, newleg)
      }
      legwrite=legwrite[!is.na(legwrite)]
      legwrite=rev(legwrite)
      
      for (i in 1:length(mod$parameters$mean)){
        print(paste("loop", i))
        print(mod$G)
        
        print(colors[i])
        #legwrite=c("G1", "G21")
        print(legwrite)
        newcolor<-append(newcolor, colors[i])
        
        p1<-p1+stat_function(fun=adjusteddnorm, n=1000, args=list(mm=mod$parameters$mean[i], sd=mod$sd[i], prop=mod$parameters$pro[i]), aes(fill=!!colors[i], color=!!colors[i]), geom="area", alpha = .2)+
          ggtitle(paste0(basename(rcmappath), "\n labels ", startlabelid, "-", endlabelid))+
          xlab("Distance (bp) between Labels of Interest")+scale_color_discrete(guide=FALSE)+ylab("Number of Molecules")+
          theme(plot.title = element_text(hjust=0.5))
    } 
    
    print("do plots")
    #print(p1)
    
    for (j in 2:length(ggplot_build(p1)$data)){
      colori<-ggplot_build(p1)$data[[j]]$colour
      print(j)
      print(unique(colori))
      newcolor<-unique(append(newcolor, colori))
    }
    newcolor<-newcolor[!is.na(newcolor)]
    print(p1+scale_fill_manual("Components\nInfo", values=newcolor, labels=legwrite)+scale_color_manual(values=newcolor, labels=legwrite, guide=FALSE))
    #return(p1)
    }
  
    while (TRUE %in% (dis$Freq < round(0.15*sum(dis$Freq),0))) {
      ###If one of the cluster is with <15% of mol
      print(paste(k, "clusters were found and needs to be reduced"))
      Gtoberemoved<-which((dis$Freq < 0.15*sum(dis$Freq)) ==TRUE)
      kremove<-dis[dis$Freq < 0.15*sum(dis$Freq),]$Var1
      k<-k-length(kremove)
      if (TRUE %in% ((which(mod$parameters$variance$sigmasq==max(mod$parameters$variance$sigmasq))) == Gtoberemoved)){
        ###If the one cluster having <15% of mol is the cluster with the largest variance
        ###Then remove data of that cluster
        filtereddata<-completeframe[!mod$classification %in% Gtoberemoved,]
      } else {
        filtereddata<-completeframe
      }
    
      #filtereddata<-completeframe
      print(paste(nrow(completeframe)-nrow(filtereddata), "data is removed"))
      mod <- Mclust(filtereddata$distance, modelNames="V", G=k)
      j<-summary(mod)
      summary(mod$BIC)
      mclust1Dplot(data = filtereddata$distance, what = "density", parameters=mod$parameters, z=mod$z)
      abline(v=mod$parameters$mean, lty=3)
      mclust1Dplot(data = filtereddata$distance, what = "classification", parameters=mod$parameters, z=mod$z)
      #mclust1Dplot(data = completeframe$distance, what = "uncertainty", parameters=mod$parameters, z=mod$z)  ###uncertainty is how uncertain it is to be asscoiated with the group classified. 
      plot.mclustBIC(mod$BIC)
      mod$sd<-sqrt(mod$parameters$variance$sigmasq)
      print(nrow(filtereddata))
      print(paste(completeframe$refdistance[1], "ref distance"))
      print(mod$parameters$mean)
      print(mod$parameters$pro)
      print(mod$sd)
      print(mod$parameters$variance$sigmasq)
      hist(mod$uncertainty)
      dis<-data.frame(table(mod$classification))
      completeframe<-filtereddata
      Plotting_hist_and_Gau(mod)
    } 
    ###Plotting begins here

   
  if ((TRUE %in% (diff(sort(mod$parameters$mean, decreasing=FALSE)) <max(sqrt(mod$parameters$variance$sigmasq)))) & (TRUE %in% (dis$Freq < 0.10*sum(dis$Freq)))){
    ###If the difference between means is smaller than the max of the std dev & that one of the clusters is less than 10% of molecules. 
    
    print(paste(k, "clusters were found and needs to be reduced"))
    k<-k-1
    filtereddata<-completeframe
    mod <- Mclust(filtereddata$distance, modelNames="V", G=k)
    j<-summary(mod)
    summary(mod$BIC)
    mclust1Dplot(data = filtereddata$distance, what = "density", parameters=mod$parameters, z=mod$z)
    abline(v=mod$parameters$mean, lty=3)
    mclust1Dplot(data = filtereddata$distance, what = "classification", parameters=mod$parameters, z=mod$z)
    
    plot.mclustBIC(mod$BIC)
    mod$sd<-sqrt(mod$parameters$variance$sigmasq)
    print(nrow(filtereddata))
    print(paste(completeframe$refdistance[1], "ref distance"))
    print(mod$parameters$mean)
    print(mod$parameters$pro)
    print(mod$sd)
    print(mod$parameters$variance$sigmasq)
    hist(mod$uncertainty)
    dis<-data.frame(table(mod$classification))
    completeframe<-filtereddata
    Plotting_hist_and_Gau(mod)
  } else {
    print("I am inside 1")
    mod$sd<-sqrt(mod$parameters$variance$sigmasq)
    ###Plotting begins here
    Plotting_hist_and_Gau(mod)
  }
    result<-c(r$CMapId[1], startlabelid, completeframe$refdistance[1], length(mod$data), k, mod$parameters$mean, mod$sd, mod$parameters$pro)
    return(result)
}
}

