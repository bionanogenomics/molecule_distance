# simplified molecule distance script 
# calculate the distance between two labels in the specified input files
# output a plot of the distances

#' @param rcmappath : anchor/reference cmap
#' @param xmap: xmap with alignment information
#' @param qcmappath: query/molecule cmap
#' @param startlabelid: ID of left reference label of interest
#' @param endlabelid: ID of right reference label of interest
#' @param out_handle: prefix for the output files
#' @return A plot of distances between labels IDs (out_handle.png)
#' #' @return A plot of distances between labels IDs (out_handle.csv)


#library(mixtools)
#library(mclust)
#library(prabclus)
library(ggplot2)
#source("C://Users//jburke//OneDrive - Bionano Genomics//Documents//services_lab_projects//radboud//scripts//get_corresponding_molecule_coordinate.R")
#source("C://Users//jburke//OneDrive - Bionano Genomics//Documents//services_lab_projects//radboud//scripts//get_header.R")

# command line args
cmd_args <- commandArgs(trailingOnly=TRUE)
rcmappath <- cmd_args[1]
xmap <- cmd_args[2]
qcmappath <- cmd_args[3]
startlabelid <- cmd_args[4]
endlabelid <- cmd_args[5]
out_handle <- cmd_args[6]

##
# helper functions
##
correspondingmolcoor<-function(rownum, reflabelnum, withalignments_1){
  permolalignment<-withalignments_1[rownum,]
  eachlabel<-strsplit(permolalignment$Alignment, split="[)(]")
  wantedqrylabels<-eachlabel
  
  getperwantedqrylabel<-function(num, reflabelnum){
    perwantedqrylabels<-wantedqrylabels[[num]]
    the_align<-perwantedqrylabels[grep(reflabelnum, perwantedqrylabels)]
    the_align<-unlist(strsplit(the_align, split=". "), recursive=F)
    peralign<-strsplit(the_align, split=",")
    findcorrespond<-data.frame(matrix(unlist(peralign), length(peralign), byrow =T))
    good<-which(findcorrespond$X1 == reflabelnum)
    if(length(good)>0){
      the_align<-the_align[good]
      qrylabel<-do.call(rbind, strsplit(the_align, split=","))[,2]
      
      return(qrylabel)
    } else {
      return("NA")
    }
  }
  qrylabels<-lapply(seq(1, length(wantedqrylabels), 1), getperwantedqrylabel, reflabelnum)
  return(qrylabels)
}

getqrycoor<-function(rcmap, xmap, labelstart){

  r<-rcmap
  x<-xmap
  
  refstart<-r[r$SiteID == labelstart,]
  refstartpos<-refstart$Position
  refstartmolcount<-refstart$Occurrence
  withalignments<-x[grepl(labelstart, x$Alignment),]
  print("var initialized")
  
  print(str(refstart))
  print(str(refstartpos))

  num_row<-as.integer(nrow(withalignments))  
  qrycoor<-lapply(seq(1, num_row, 1), correspondingmolcoor, labelstart, withalignments_1=withalignments)
  print("qrycoor")
  
  result<-data.frame(withalignments$QryContigID, do.call(rbind, qrycoor))
  print("result data frame initiated")
  print(str(result))
  
  result$refPos<-refstartpos
  print("refPos update")
  print(str(result))
  
  result$Qrylabel<-unlist(result$Qrylabel)
  print("qrylabel update 1")
  print(str(result))
  
  colnames(result)<-c("QryContigID", "Qrylabel", "refpos")
  result$Qrylabel<-unlist(result$Qrylabel)
  print("qrylabel update 2")
  print(str(result))
  
  result<-result[!result$Qrylabel=="NA",]
  print(str(result))
  
  return(result)
}

getcoor<-function(n, starts, qcmap_1){
  perqry<-starts[n,]
  findcoor<-qcmap_1[qcmap_1$CMapId == perqry$QryContigID,]
  findcoor<-findcoor[findcoor$SiteID == perqry$Qrylabel,]$Position
  return(findcoor)
}

getcolnames<-function(rarepath){
  rare1<-read.table(rarepath, comment.char="#", sep="\t", stringsAsFactors = F)
  nhl<-system(paste0("grep '^#' ", rarepath, "| wc -l"), intern =TRUE)
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
  return(rare1)
}

getheader<-function(rarepath){
  rare1<-read.table(rarepath, comment.char="#", sep="\t", stringsAsFactors = F)
  nhl<-system(paste0("grep '^#' ", rarepath, " | wc -l"), intern =TRUE)
  header<-readLines(rarepath, n=nhl)
  return(header)
}

waterfall_plotting <- function(distance_data, filename_prefix){
  distance_mean <- mean(distance_data$distance)
  distance_data$index <- row.names(distance_data)
  waterfall_plot<- ggplot(distance_data, aes(x=reorder(index,distance), y=distance)) +
    geom_bar(stat="identity", 
             fill="#0072B2",
             width=0.7,
             position = position_dodge(width=0.4))+
    scale_y_discrete(expand = c(0,0)) +
    labs(title = "Distance Between Specified Labels", 
         x = "Molecules", 
         y = "Distance (basepairs)") +
    theme_classic()+
    theme(axis.text.x = element_blank()) +
    theme(legend.position="none") +
    geom_hline(yintercept=distance_mean, linetype="dashed", color = "#D55E00")+
    geom_text(aes(10,distance_mean,label = "mean", vjust = -1), 
             color = "#D55E00", check_overlap=T )
    ggsave(plot = waterfall_plot, file = paste0(filename_prefix,'_barplot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
  }

summary_stats <- function(distance_data, filename_prefix){
  mean_distance<-mean(distance_data$distance)
  sd_distances<- sd(distance_data$distance)
  cv_distances <- sd_distances/mean_distance
  summary_table <- cbind(mean_distance, sd_distances,cv_distances)
  write.csv(file=paste0(filename_prefix, '_summary_stats.csv'), x=summary_table)
  print(summary_table)
}

histogram_plot <- function(distance_data, filename_prefix){
  num_bins <- sqrt(nrow(distance_data))
  hist_dens_plot <- ggplot(distance_data, aes(x = distance)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "bisque1", bins=num_bins) +
    geom_density()+
    labs(title = "Histogram of Distances with Density Estimate", 
         x = "Distance (basepairs)", 
         y = "Density") +
    theme_classic()
    ggsave(plot = hist_dens_plot, file = paste0(filename_prefix,'_histogram.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300) 
}

violin_plot <- function(distance_data, filename_prefix){
  distance_data$arbitrary_group <- 'xyz'
  violin_with_boxplot <- ggplot(distance_data, aes(y=distance, x=arbitrary_group)) + 
    geom_violin(trim=FALSE, fill="darkslategray4")+
    labs(title="Violin Plot of Distances (with Boxplot)", x = "Sample", y = "Distances (basepairs)")+
    geom_boxplot(width=0.1)+
    theme_classic()+ 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
    ggsave(plot = violin_with_boxplot, file = paste0(filename_prefix,'_violin_plot.png'), height = 5, width= 7, units = "in", device = 'png', dpi=300)
}

### "main" function
molecule_distance_main<-function(rcmappath,
                                 xmap,
                                 qcmappath,
                                 startlabelid, 
                                 endlabelid, 
                                 out_handle
                                 #additional_noise_removal=FALSE # is this needed bc I removed it; was in Joyce's script but we've always had it "off" (e.g. set to FALSE)
                                 )
{
  print("INPUT FILES:")
  print(rcmappath)
  print(xmap)
  print(qcmappath)
  qcmap<-getcolnames(qcmappath)
  qcmap$CMapId<-as.character(qcmap$`#h CMapId`)
  print(rcmappath)
  r<-getcolnames(rcmappath)
  r$CMapId<-as.character(r$`#h CMapId`)
  x<-read.table(xmap, stringsAsFactors=F)
  h<-getheader(xmap)
  h2<-h[length(h)-1]
  h<-strsplit(h2, split="\t")[[1]]
  colnames(x)<-h
  
  starts<-getqrycoor(r, x, startlabelid)
  ends<-getqrycoor(r, x, endlabelid)
  
  ### Downsample to the most 1000 molecules
  if(dim(starts)[1]>1000){
    starts<-starts[sample(nrow(starts),1000),]
    ends<-ends[sample(nrow(ends), 1000),]
  }
  
  ### get the righmost label if the mols that have multiple alignments to one ref label
  starts$Qrylabel<-unlist(lapply(starts$Qrylabel, max))
  ### get the leftmost label if the mols have multiple alignments to one ref label
  ends$Qrylabel<-unlist(lapply(ends$Qrylabel, min))
  
  is.error <- function(x) inherits(x, "try-error")
  
  ### Get Qry position
  if(nrow(starts)<10 | nrow(ends)<10){   ###Too few molecules at the start of end (Added on 04/20/2021)
    result<-c(r$CMapId[1], startlabelid, starts$refpos[1], min(nrow(starts), nrow(ends)), 0, 0, 0, 0)
    print(paste("Only <10 molecules are the label start or label end"))
    return (result)
  }
  
  starts$QryPos<-lapply(seq(1, nrow(starts), 1), getcoor, starts, qcmap_1=qcmap)
  ends$QryPos<-lapply(seq(1, nrow(ends),1), getcoor, ends, qcmap_1=qcmap)
  completeframe <- merge(starts,ends,by="QryContigID")
  print(paste0(nrow(starts)," molecules at start, ", nrow(ends), " molecules at end, ", nrow (completeframe), " spanned both labels"))
  
  ### x is left breakpoint, y is right breakpoint
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
    simple_output <- completeframe[c("QryContigID","distance")]
    write.csv(file = paste0(out_handle, "_distances.csv"), x = simple_output)
    write.csv(file = paste0(out_handle, "_complete_data.csv"), x = completeframe)
    print("QC pass")
  }
  waterfall_plotting(simple_output, out_handle)
  histogram_plot(simple_output, out_handle)
  violin_plot(simple_output, out_handle)
  summary_stats(simple_output, out_handle)
}

# function call
molecule_distance_main(rcmappath, xmap, qcmappath, startlabelid, endlabelid, out_handle)

