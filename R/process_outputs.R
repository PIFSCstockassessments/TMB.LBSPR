#' TMB.LPSBR model output
#'
#' Takes the output from run_analyses, a Monte Carlo dataset, and creates tables and graphs
#'
#' @param D TMB.LBSPR MCMC dataset
#' @param TYPE ype of data available ("Survey only", "Catch only", "Both).
#' @param SHOW.LC Option to show minimum size (Lc) management results
#' @param outdir Location to store outputs to file. Defaults to TMB.LBSPR directory in the system's Home.
#'
#' @import ggplot2
#' @import scales
#' @import gridExtra
#' @import grid
#' @import MASS
#' @import Hmisc
#' @import openxlsx
#' @importFrom reshape2 melt dcast
#' @importFrom data.table data.table
#'
process_outputs <- function(D, TYPE=c("Both", "Survey only", "Catch only"), SHOW.LC=TRUE,
                            outdir=file.path(home=Sys.getenv("HOME"), "TMB.LBSPR")){


  TYPE <- match.arg(TYPE)

  options(na.rm=TRUE)
  D <- data.frame(D)

  #TYPE<-"Both"   # COMMENT OUT, JUST FOR TESTING
  #SHOW.LC<-T     # SAME
  #D <- Final     # SAME

#=====================FILTERS=============================================
D <- D[D$Bio.survey>0,]
D <- D[D$Linf>0,]
D <- D[D$F30>-9999,]
D <- D[D$Fmort>0.005,]

#if(TYPE!="Survey only"){D <- D[C30_Catch>0]}
#if(TYPE!="Catch only") {D <- D[C30_Survey>0]}
D[D$Lc30<=50|D$Lc30==-9999,]$Lc30 <- NA
#D <- data.frame(D); D <- data.table(D)
#D[D$Lc30==-9999,]$Lc30 <- NA
#D <- data.frame(D)
D <- data.frame(D)

# Calculate selectiviy medians
MEDIAN  <- prettyNum(signif(sapply(D[,14:15],median,na.rm=T),3))
SD      <- prettyNum(signif(sapply(D[,14:15],sd,na.rm=T),3))
Summary.select <- cbind(MEDIAN,SD)

D <- D[,c("Linf","K","Lmat50","Amat","Amax","M","F30","Lc30","Catch","Bio.catch","C30.catch","Bio.survey","C30.survey","Fmort","SPR","F_Fmsy")]

#=====================TRANSFORMATIONS================================
for(j in 9:13){D[,j]<-D[,j]/1000} #kg to 1000 kgs

#====================GRAPH METADATA INFORMATION===========================
Min  <- vector(length=ncol(D))
Max  <- vector(length=ncol(D))
Units <- c("mm","yr^-1","mm","yr","yr","yr^-1","yr^-1","mm","kg","kg","kg","kg","kg","yr^-1","","")

Titles <- list(expression(paste(italic("L")["inf"]~(mm))),
               expression(paste(italic("K")~(yr^-1))),
               expression("Length at maturity"~(mm)),
               expression("Age at maturity"~(mm)),
               expression("Longevity"~(yr)),
               expression(paste(italic("M")~(yr^-1))),
               expression(paste(italic("F")["30"]~(yr^-1))),
               expression(paste(italic("Lc")["30"]~(mm))),
               expression("Total catch"~("1000kg")),
               expression(paste(italic("B")~"from catch"~("1000kg"))),
               expression(paste(italic("C")["30"]~"from catch"~("1000kg"))),
               expression(paste(italic("B")~"from divers"~("1000kg"))),
               expression(paste(italic("C")["30"]~"from divers"~("1000kg"))),
               expression(paste(italic("F")~(yr^-1))),
               expression(paste(italic("SPR"))),
               expression(paste(italic("F")/italic("F")["30"]))
)

create.expr <- function(num, unit ){

  if(unit!="")
    text <- paste("expression(",'"', num,'"', "~", unit, ")", sep="")
  else
    text <- paste("expression(",'"', num,'"',")", sep="")
  express <- eval(parse(text=text))
  return (express)
}

temp <- D; temp[temp==-9999]<-NA
for (i in 1:ncol(temp)){ # Find min/max of each columns
  Min[i]   <- min(temp[,i],na.rm=T)
  Max[i]   <- max(temp[,i],na.rm=T)
}

if(TYPE=="Both"){
  Min[10] <- Min[12]  <- min(temp[,10],temp[,12])
  Max[10] <- Max[12]  <- max(mean(temp[,10])+sd(temp[,10]), mean(temp[,12])+sd(temp[,12])*6)
  Min[11] <- Min[13]  <- min(temp[,11],temp[,13])
  Max[11] <- Max[13]  <- max(mean(temp[,11])+sd(temp[,11]),mean(temp[,13])+sd(temp[,13])*6)

} else if(TYPE=="Survey only"){
  Min[12]<- min(temp[,12])
  Max[12]<- mean(temp[,12])+sd(temp[,12])*6

  Min[13]<- min(temp[,13])
  Max[13]<- mean(temp[,13])+sd(temp[,13])*6

} else if(TYPE=="Catch only"){
  Min[10]<- min(temp[,10])
  Max[10]<- mean(temp[,10])+sd(temp[,10])*2

  Min[11]<- min(temp[,11])
  Max[11]<- mean(temp[,11])+sd(temp[,11])*2
}

Min[15] <- 0; Max[15]<-1

#======================GRAPH FORMATING=====================================
#colnames(Meta) <- c("Par","Name","Unit1","Unit2","Min","Max")
theme <- theme(axis.line=element_line(colour="black"),
               panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               panel.border=element_blank(),
               panel.background=element_blank(),
               axis.text=element_text(color="black"),
               axis.title = element_text(size=10)
)

#=====================CREATE GRAPHS========================================
glist <- vector(mode="list",length=17)
for(i in 1:16){

  if(SHOW.LC==F&i==8){next}

  bw       <- diff(range(D[,i],na.rm=T)) / (  diff(range(D[,i],na.rm=T)) / (2 * IQR(D[,i],na.rm=T) / length(D[,i])^(1/3))   )
  maxCount <- max(hist(D[,i],breaks=seq(min(D[,i],na.rm=T),max(D[,i],na.rm=T)+bw,by=bw),plot=F)$counts)
  max.X    <- max(hist(D[,i],breaks=seq(min(D[,i],na.rm=T),max(D[,i],na.rm=T)+bw,by=bw),plot=F)$breaks)
  min.X    <- min(hist(D[,i],breaks=seq(min(D[,i],na.rm=T),max(D[,i],na.rm=T)+bw,by=bw),plot=F)$breaks)
  range.X  <- max.X - min.X
  TwenPerc <- 0.25*range.X

  aMedian  <- median(D[,i],na.rm=T)
  if(i==8){ # This section is to include NAs when calculating the Lc30 median
    D8             <- D[,8]
    D8[is.na(D8)]  <- 0
    aMedian        <- median(D8)
   }

  if(i==15) {
    bw <- 0.02
    maxCount <- max(hist(D[,i],breaks=seq(min(D[,i],na.rm=T),max(D[,i],na.rm=T)+bw,by=bw),plot=F)$counts)
    max.X    <- max(hist(D[,i],breaks=seq(min(D[,i],na.rm=T),max(D[,i],na.rm=T)+bw,by=bw),plot=F)$breaks)
    min.X    <- min(hist(D[,i],breaks=seq(min(D[,i],na.rm=T),max(D[,i],na.rm=T)+bw,by=bw),plot=F)$breaks)
    range.X  <- max.X - min.X
    }

  MedianLabel<-""

  Median.pos <- aMedian+bw*6
  if(i==8){Median.pos <- aMedian-TwenPerc}

  if(i>=9&i<=13){MedianLabel<- prettyNum(signif(aMedian*1000,3),big.mark=",")}else{MedianLabel<- prettyNum(signif(aMedian,2),big.mark=",") }

  anExpr <- create.expr(MedianLabel,Units[i])

  variable <- colnames(D)[i]
  aGraph <- ggplot(data=D,aes_string(x=variable))+geom_histogram(binwidth=bw,col="cadetblue3",fill="cadetblue3",lwd=0.2)+
    scale_x_continuous(name=Titles[[i]],labels=comma)+
    coord_cartesian(xlim=c(as.numeric(Min[i]),as.numeric(Max[i])))+
    scale_y_continuous(expand=c(0,0))+
    geom_vline(xintercept=aMedian,col="red",linetype="dashed",size=1)+
    annotate("text",hjust=0,x=Median.pos,y=maxCount*0.8,label=as.character(anExpr),size=2.5,parse=T)+
    theme_bw()+theme


  if(i==15){
    aGraph <- aGraph + annotate("segment",x=0.3,xend=0.3,y=0,yend=maxCount/20,col="black",size=1)
  }

  glist[[i]] <- aGraph

}

#========================ARRANGE GRAPH PAGES==================================
LH.page     <- arrangeGrob(ncol=2, glist[[1]],glist[[2]],glist[[3]],glist[[4]],glist[[6]],glist[[5]])

if(SHOW.LC)
  Status.page <- arrangeGrob(ncol=2, glist[[14]],glist[[7]],glist[[16]],glist[[8]],glist[[15]])

if(!SHOW.LC)
  Status.page <- arrangeGrob(ncol=2, glist[[14]],glist[[7]],glist[[16]],glist[[15]])

if(TYPE=="Both")
  #C30.page    <- arrangeGrob(ncol=2, glist[[9]],glist[[11]],glist[[10]],glist[[13]],glist[[12]])
  C30.page    <- arrangeGrob(ncol=2, glist[[11]],glist[[10]],glist[[13]],glist[[12]],glist[[9]])

if(TYPE=="Survey only")
  C30.page    <- arrangeGrob(ncol=2, glist[[13]],glist[[12]]  )

if(TYPE=="Catch only")
  C30.page    <- arrangeGrob(ncol=2, glist[[11]], glist[[10]],glist[[9]]  )

#grid.draw(LH.page)
#grid.draw(Status.page)
#grid.draw(C30.page)


#==========================CREATE CUMULATIVE C30 PLOT==========================================================

cumul            <- NULL
cumul$prob       <- seq(0,0.95,by=0.01)
cumul$C30.catch  <- quantile(D$C30.catch, probs=seq(0,0.95,by=0.01))
cumul$C30.survey <- quantile(D$C30.survey, probs=seq(0,0.95,by=0.01))
cumul            <- data.frame(cumul)

C30.cdf.graph <- ggplot(data=cumul,aes_(y=quote(prob)))

if(TYPE=="Both"){
  C30.cdf.graph <- C30.cdf.graph+geom_line(aes(x=cumul$C30.catch),col="orange",size=0.8,linetype="longdash")+
    geom_line(aes(x=cumul$C30.survey),col="blue",size=0.8,linetype="dotted")+
    annotate("segment",x=median(D$C30.catch),xend=median(D$C30.catch),y=0,yend=0.05,col="orange",size=1)+
    annotate("segment",x=median(D$C30.survey),xend=median(D$C30.survey),y=0,yend=0.05,col="blue",size=1)


} else if(TYPE=="Catch only") {
  C30.cdf.graph <- C30.cdf.graph+geom_line(aes(x=cumul$C30.catch),col="orange",size=0.8,linetype="longdash")+
    annotate("segment",x=median(D$C30.catch),xend=median(D$C30.catch),y=0,yend=0.05,col="orange",size=1)
} else {
  C30.cdf.graph <- C30.cdf.graph+geom_line(aes(x=cumul$C30.survey),col="blue",size=0.8,linetype="dotted")+
    annotate("segment",x=median(D$C30.survey),xend=median(D$C30.survey),y=0,yend=0.05,col="blue",size=1)
}

C30.cdf.graph<- C30.cdf.graph+scale_x_continuous(name=expression(paste(italic("C")[30]~("1000kg"))),expand=c(0,0))+
  scale_y_continuous(name="Overfishing prob.",expand=c(0,0),limits=c(0,.95))+
  theme_bw()+theme


#===============================SAVE GRAPHS TO FILE=====================================
if(TYPE=="Both"){height=16}else if(TYPE=="Catch only"){height=11}else if(TYPE=="Survey only"){height=5.5}

filename <- file.path(outdir,"1_LH_page.tiff")
tiff(filename=filename, type="cairo", units="cm", compression = "lzw",
     width=14,
     height=17,
     res=400)
grid.draw(LH.page)
dev.off()

if(SHOW.LC){height=14}else if(!SHOW.LC){height=9.33}
filename <- file.path(outdir,"2_Status_page.tiff")
tiff(filename=filename, type="cairo",units="cm", compression = "lzw",
     width=12.5,
     height=height,
     res=400)
grid.draw(Status.page)
dev.off()

if(TYPE=="Both"){height=17}else if(TYPE=="Catch only"){height=11.33}else if(TYPE=="Survey only"){height=5.666}
filename <- file.path(outdir,"3_C30_page.tiff")
tiff(filename=filename, type="cairo",units="cm", compression = "lzw",
     width=14,
     height=height,
     res=400)
grid.draw(C30.page)
dev.off()

filename <- file.path(outdir,"4_C30_cdf.tiff")
tiff(filename=filename, type="cairo",units="cm", compression = "lzw",
     width=14,
     height=6,
     res=400)
print(C30.cdf.graph)
dev.off()


#=========GENERATE A PROBABILITY OF OVERFISHING TABLE FOR VARIOUS C30s=====
C30.table <- NULL
C30.table$Overfishing_prob  <- seq(0.05,0.5,by=0.01)
C30.table$C30.catch         <- NULL
C30.table$C30.survey        <- NULL

if(TYPE=="Both"){
  C30.table$C30.catch  <- quantile(D$C30.catch, probs=seq(0.05,0.5,by=0.01))
  C30.table$C30.survey <- quantile(D$C30.survey, probs=seq(0.05,0.5,by=0.01))
  C30.table$C30.survey <- round(C30.table$C30.survey,3)
  C30.table$C30.catch  <- round(C30.table$C30.catch,3)
  }
if(TYPE=="Survey only"){
  C30.table$C30.survey <- quantile(D$C30.survey, probs=seq(0.05,0.5,by=0.01))
  C30.table$C30.survey <- round(C30.table$C30.survey,3)
  }

if(TYPE=="Catch only"){
  C30.table$C30.catch  <- quantile(D$C30.catch, probs=seq(0.05,0.5,by=0.01))
  C30.table$C30.catch <- round(C30.table$C30.catch,3)
}

C30.table <- data.frame(C30.table)
C30.table1 <- C30.table[C30.table$Overfishing_prob<0.28,]
C30.table2 <- C30.table[C30.table$Overfishing_prob>=0.28,]
C30.table3 <- cbind(C30.table1,C30.table2)

C30.table <- format(C30.table,digits=3)
C30.table$Overfishing_prob <- format(C30.table$Overfishing_prob,digits=2)

wb <- createWorkbook()
addWorksheet(wb,"C30")
writeData(wb,sheet="C30",C30.table3)
style <- createStyle(halign="center",fontName="Times New Roman")
addStyle(wb,sheet="C30",style,cols=c(1:6),rows=1:50,gridExpand=T)

#===GENERATE A PROBABILITY OF OVERFISHING TABLE FOR VARIOUS LC=====
LC30.table <- NULL
LC30.table$prob <- seq(0.5,0.95,by=0.01)
LC30.table <- data.frame(LC30.table)

D_Lc30                <- D$Lc30
D_Lc30[is.na(D_Lc30)] <- 0

LC30 <- NULL
LC30 <- quantile(D_Lc30,probs=seq(0.5,0.95,by=0.01),na.rm=T)

LC30.table$LC30 <- round(LC30,0)
LC30.table$prob <- 1-LC30.table$prob
LC30.table      <- LC30.table[order(LC30.table$prob),]

LC30.table       <- data.frame(LC30.table)
LC30.table1      <- LC30.table[LC30.table$prob<0.28,]
LC30.table1$prob <- format(LC30.table1$prob,digits=2)
LC30.table2      <- LC30.table[LC30.table$prob>=0.28,]
LC30.table2$prob <- format(LC30.table2$prob,digits=2)
LC30.table3      <- cbind(LC30.table1,LC30.table2)

LC30.table3      <- format(LC30.table3,digits=0)
LC30.table3$prob <- format(LC30.table3$prob,digits=2)


addWorksheet(wb,"LC30")
writeData(wb,sheet="LC30",LC30.table3)
addStyle(wb,sheet="LC30",style,cols=c(1:6),rows=1:50,gridExpand=T)

filename <- file.path(outdir,"MANAGEMENT.xlsx")
saveWorkbook(wb,filename,overwrite=T)


#=====SUMMARY TABLE FOR REPORT===========================
D[is.na(D$Lc30),]$Lc30 <- 0
MEDIAN  <- prettyNum(signif(sapply(D,median,na.rm=T),3))
SD      <- prettyNum(signif(sapply(D,sd,na.rm=T),3))
Summary <- cbind(MEDIAN,SD)

# Calculate iteration below SPR 0.3
D <- data.table(D)
num_iter <- round(nrow(D[D$SPR<0.3,])/nrow(D),3)

Summary <- rbind(Summary,num_iter,Summary.select)

filename <- file.path(outdir,"Summary.csv")
write.csv(Summary,file=filename)

#======================CALCULATE THE PERCENT SPR ITERATIONS BELOW 30%==================


} # End of function









