###
PERMANOVAandANOSIM=function(TT,Meta){
  library("vegan")
  library("reshape2")
  tt=t(TT)
  Meta_CE1Only=Meta[which(sapply(strsplit(rownames(Meta), "H" , fixed = TRUE), "[", 1)=="CE1"),]
  Meta_CE1Only=rbind(Meta_CE1Only,
                     Meta[which(rownames(Meta)=="CE1XenicR1"),],
                     Meta[which(rownames(Meta)=="CE1XenicR2"),])
  
  META=Meta_CE1Only[rownames(tt),]
  all(rownames(tt) == rownames(META))
  
  permanova=adonis2(tt ~ Hratio, data = META, permutations = 999, method="bray", by = NULL)  
  RPERMANOVA=permanova$R2[1]
  pPERMANOVA=permanova$`Pr(>F)`[1]
  if (pPERMANOVA<=0.001){
    pPERMANOVA="<0.001"
  }else{
    pPERMANOVA=pPERMANOVA
  }
  PermanovaTable=permanova
  ####ANOSIM###
    MM=matrix(0,ncol(TT),ncol(Meta))
    MM=list()
    NM=as.numeric()
    for (i in 1:ncol(TT)){
      M=Meta
      MM[[i]]=M[which(rownames(M)==colnames(TT)[i]),][5]
      NM[i]=rownames(M[which(rownames(M)==colnames(TT)[i]),][5])
    }
    Hratio=melt(MM)
    rownames(Hratio)=NM
    BC_dist = avgdist(t(TT),2000,  dmethod = "bray")
    dune.ano = with(TT, anosim(BC_dist, Hratio[,1]))
    RANOSIM=dune.ano$statistic
    pANOSIM=dune.ano$signif
    if (pANOSIM<=0.001){
      pANOSIM="<0.001"
    }else{
      pANOSIM=pANOSIM
    }
  
  PERMANOVAandANOSIM = vector('expression',4)
  PERMANOVAandANOSIM[1] = substitute(expression(italic(ANOSIM-R) == MYVALUE), 
                         list(MYVALUE = format(RANOSIM, digits = 2)))[2]
  PERMANOVAandANOSIM[2] = substitute(expression(italic(ANOSIM-p) == MYVALUE), 
                                     list(MYVALUE = format(pANOSIM, digits = 2)))[2]
  PERMANOVAandANOSIM[3] = substitute(expression(italic(PERMANOVA-R^2) == MYVALUE), 
                                     list(MYVALUE = format(RPERMANOVA, digits = 2)))[2]
  PERMANOVAandANOSIM[4] = substitute(expression(italic(PERMANOVA-p) == MYpVALUE), 
                         list(MYpVALUE = format(pPERMANOVA, digits = 2)))[2]

  legend('bottomleft', legend = PERMANOVAandANOSIM, bty = 'n',cex=0.6)
 return(PermanovaTable)
}

###
OrdinationEachDay=function(MyData,day){
  DDay=filter(MyData, MyData$Day==day)
  plot(DDay[,1],DDay[,2],pch=1, col="white", 
       xlim=c(min(MyData[,1]),max(MyData[,1])),
       ylim=c(min(MyData[,2]),max(MyData[,2])),
       type="b",xaxt="n",yaxt="n", 
       xlab="", ylab="")
  DD=filter(DDay, DDay$H==0)
  points(DD[,1],DD[,2],pch=1, col="black")
  DD=filter(DDay, DDay$H==5)
  points(DD[,1],DD[,2],pch=2, col="black")
  DD=filter(DDay, DDay$H==25)
  points(DD[,1],DD[,2],pch=16, col="gray")
  DD=filter(DDay, DDay$H==50)
  points(DD[,1],DD[,2],pch=17, col="gray")
  DD=filter(DDay, DDay$H==75)
  points(DD[,1],DD[,2],pch=16, col="black")
  DD=filter(DDay, DDay$H==100)
  points(DD[,1],DD[,2],pch=17, col="black")
  DD=filter(MyData, rownames(MyData)=="CE1XenicR1")
  points(DD[,1],DD[,2],pch=8, col="black")
  DD=filter(MyData, rownames(MyData)=="CE1XenicR2")
  points(DD[,1],DD[,2],pch=8, col="black")
  return(DDay)
}

Tukeytest=function(model,a){
  TUKEY=TukeyHSD(model, conf.level=.95) 
  #plot(TukeyHSD(model, conf.level=.95), las = 2)
  # I want to write the letter over each box. Over is how high I want to write it.
  #Add the labels
  Tukey.levels = TUKEY[[1]][,4]
  Tukey.labels = data.frame(multcompLetters(Tukey.levels)['Letters'])
  Tukey=Tukey.labels[order(as.numeric(rownames(Tukey.labels))),]
  
  return(Tukey)
}

ANOVAmodelandPlot=function(a,model,Tukey,x,y,ylim){
  p=summary(model)[[1]][1,5]
  if (p<=0.05){
    Tukey=Tukey
  }else{
    Tukey=rep("",6)
  }
 
  t=(ylim[2]-ylim[1])*0.075
  text(c(1:6), a$stats[nrow(a$stats),]+t,Tukey, cex=0.75)
  
  if (p<=0.001){
    p="<0.001"
  }else{
    p=p
  }
  ANOVAp = vector('expression',2)
  ANOVAp[1] = substitute(expression(italic(p) == MYVALUE), 
                         list(MYVALUE = format(p, digits = 2)))[2]
  text(x,y, label = ANOVAp[1], cex=0.7)
}

DiveristyPlot=function(MyData, D, col,R, B, ylim, Xaxis,LABLE){
  library(multcompView)
  library(vegan)
  library(ggpubr)
  library(phyloseq)
  library(dplyr)
  day=c(0,2,4,6,8,10,12)
  daytext=c("Day0","Day2","Day4","Day6","Day8","Day10","Day12")
  par(mar = c(0, 1.4, 0, 0.3))
  ##Day2
  DDay=filter(MyData, Day==day[2])
  a=boxplot(DDay[,D]~as.numeric(DDay$H),col=col,
            xaxt="n",yaxt="n",ylim=c(min(DDay[,D]),max(DDay[,D]*1.1)))
  a
  stripchart(DDay[,D]~as.numeric(DDay$H), col=281, pch=1, vertical = TRUE,add = TRUE)     
  #text(x=1.6, 31, label = daytext[2], cex = 1)
  text(x=0.7, y = max(DDay[,D]*1.1)*0.97, label = LABEL[1], cex = 1)
  axis(2, tck = -0.04,at = seq(R[1],R[2], by = B), 
       label = rep("",length(seq(R[1],R[2], by = B))))
  axis(2, las=1, cex.axis=0.7, line = -0.6, lwd = 0)

  #fit one-way ANOVA model and do Tukey test
  model=aov(DDay[,D]~DDay$H, data=DDay)
  Tukey=Tukeytest(model,a)
  ANOVAmodelandPlot(a,model,Tukey,5.6,max(DDay[,D]*1.1)*0.97,c(min(DDay[,D]),max(DDay[,D])))
  if (Xaxis==TRUE) {
    axis(1, at = c(1,2,3,4,5,6), 
         label = rep("", 6), 
         tck = -0.04)
    axis(1, at = c(1,2,3,4,5,6), 
         label = c("0%","5%","25%","50%","75%","100%"), 
         cex.axis=0.6,line = -0.7, lwd = 0)
  }
  
  ###Day4
  par(mar = c(0, 1.4, 0, 0)) # make the plots be closer together
  
  DDay=filter(MyData, Day==day[3])
  a=boxplot(DDay[,D]~as.numeric(DDay$H),col=col,
            xaxt="n",yaxt="n",ylim=ylim)
  a
  stripchart(DDay[,D]~as.numeric(DDay$H), col=281, pch=1, vertical = TRUE,add = TRUE)     
  #text(x=1.6, y =40.8, label = daytext[3], cex = 1)
  text(x=0.7, y = ylim[2]*0.97, label = LABEL[2], cex = 1)
  axis(2, tck = -0.04, at = seq(R[1],R[2], by = B), 
       label = rep("",length(seq(R[1],R[2], by = B))))
  axis(2, las=1, cex.axis=0.7, line = -0.6, lwd = 0,
       at = seq(R[1],R[2], by = B), 
       label = seq(R[1],R[2], by = B))

  #fit one-way ANOVA model and do Tukey test
  model=aov(DDay[,D]~DDay$H, data=DDay)
  Tukey=Tukeytest(model,a)
  ANOVAmodelandPlot(a,model,Tukey,5.6,ylim[2]*0.97,ylim)
  if (Xaxis==TRUE) {
    axis(1, at = c(1,2,3,4,5,6), 
         label = rep("", 6), 
         tck = -0.04)
    axis(1, at = c(1,2,3,4,5,6), 
         label = c("0%","5%","25%","50%","75%","100%"), 
         cex.axis=0.6,line = -0.7, lwd = 0)
  }
  
  ###Day6-12

  par(mar = c(0, 0, 0, 0)) # make the plots be closer together
  for (i in 4:7){
    DDay=filter(MyData, MyData$Day==day[i])
    a=boxplot(DDay[,D]~as.numeric(DDay$H), col=col,
              xaxt="n",yaxt="n",ylim=ylim)
    a
    stripchart(DDay[,D]~as.numeric(DDay$H), col=281, pch=1, vertical = TRUE, add = TRUE)    
    #text(x=1.6, y = 40.8, label = daytext[i], cex = 1)
    text(x=0.7, y = ylim[2]*0.97, label = LABEL[i-1], cex = 1)
    #fit one-way ANOVA model and do Tukey test
    model=aov(DDay[,D]~DDay$H, data=DDay)
    Tukey=Tukeytest(model,a)
    ANOVAmodelandPlot(a,model,Tukey,5.6,ylim[2]*0.97,ylim)  
    if (Xaxis==TRUE) {
      axis(1, at = c(1,2,3,4,5,6), 
           label = rep("", 6), 
           tck = -0.04)
      axis(1, at = c(1,2,3,4,5,6), 
           label = c("0%","5%","25%","50%","75%","100%"), 
           cex.axis=0.6,line = -0.7, lwd = 0)
    }
  }
}
