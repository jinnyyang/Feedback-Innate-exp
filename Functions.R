###
PERMANOVAandANOSIM=function(RawTable,Filtered_Meta,day){
  Meta_Day=filter(Filtered_Meta,Filtered_Meta$Day==day)
  RawTable_Day=matrix(0,nrow(RawTable),length(rownames(Meta_Day)))
  for (i in 1:length(rownames(Meta_Day))){
    RawTable_Day[,i]=RawTable[,which(colnames(RawTable)==rownames(Meta_Day)[i])]
  }
  rownames(RawTable_Day)=rownames(RawTable)
  colnames(RawTable_Day)=rownames(Meta_Day)
  
  BC_dist = avgdist(t(RawTable_Day), 2000,  dmethod = "bray") ## calculate Bray-Curtis distance matrix
  
  permanova=adonis2(BC_dist~Hratio, data= Meta_Day, permutations = 999, method="bray")
  RPERMANOVA=permanova$R2[1]
  pPERMANOVA=permanova$`Pr(>F)`[1]
  if (pPERMANOVA<=0.001){
    pPERMANOVA="<0.001"
  }else{
    pPERMANOVA=pPERMANOVA
  }
  PermanovaTable=permanova

  ####ANOSIM###
  dune.ano = with(Meta_Day, anosim(BC_dist, Hratio))
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
OrdinationEachDay=function(MyData,day,COL,PCH){
  DDay=filter(MyData, MyData$Day==day)
  plot(DDay[,1],DDay[,2],pch=1, col="white", 
       xlim=c(min(MyData[,1]),max(MyData[,1])),
       ylim=c(min(MyData[,2]),max(MyData[,2])),
       type="b",xaxt="n",yaxt="n", 
       xlab="", ylab="")
  DD=filter(DDay, DDay$H==0)
  points(DD[,1],DD[,2],pch=PCH[6], col=COL[6])
  DD=filter(DDay, DDay$H==5)
  points(DD[,1],DD[,2],pch=PCH[5], col=COL[5])
  DD=filter(DDay, DDay$H==25)
  points(DD[,1],DD[,2],pch=PCH[4], col=COL[4])
  DD=filter(DDay, DDay$H==50)
  points(DD[,1],DD[,2],pch=PCH[3], col=COL[3])
  DD=filter(DDay, DDay$H==75)
  points(DD[,1],DD[,2],pch=PCH[2], col=COL[2])
  DD=filter(DDay, DDay$H==100)
  points(DD[,1],DD[,2],pch=PCH[1], col=COL[1])
  DD=filter(MyData, rownames(MyData)=="CE1XenicR1")
  points(DD[,1],DD[,2],pch=PCH[7], col=COL[7])
  DD=filter(MyData, rownames(MyData)=="CE1XenicR2")
  points(DD[,1],DD[,2],pch=PCH[7], col=COL[7])
  return(DDay)
}
###
CorrelationCoefficientEachH=function(Day,CE1_Only){
  all(rownames(Day)==rownames(CE1_Only))
  ###Prepare relative abundance data
  Day=sapply(strsplit(sapply(strsplit(colnames(CE1_Only), "D" , fixed = TRUE), "[", 2), "R" , fixed = TRUE), "[", 1)
  D_CE1_Only=CE1_Only[,which(Day==12)]
  CE1_Only=CE1_Only[rowSums(CE1_Only)>0,]
  Hratio=sapply(strsplit(sapply(strsplit(colnames(D_CE1_Only), "D" , fixed = TRUE), "[", 1), "H" , fixed = TRUE), "[", 2)
  SUMReadEachHratio=data.frame(
    rowSums(D_CE1_Only[,which(Hratio==0)])/5,
    rowSums(D_CE1_Only[,which(Hratio==5)])/5,
    rowSums(D_CE1_Only[,which(Hratio==25)])/5,
    rowSums(D_CE1_Only[,which(Hratio==50)])/5,
    rowSums(D_CE1_Only[,which(Hratio==75)])/5,
    rowSums(D_CE1_Only[,which(Hratio==100)])/5)
  DATA=data.frame(Day,SUMReadEachHratio)
  DATA=DATA[which(rowSums(DATA)>0),]
  colnames(DATA)=c("Cor","P-value","SUM","H0Mean","H5Mean","H25Mean","H50Mean","H75Mean","H100Mean")
  return(DATA)
}

###

figure5=function(MyData,D,a,b){
  plot(as.numeric(MyData$H),MyData[,D],ylim=c(a,b),xaxt="n", 
       xlim=c(-5,105))
  l=lm(MyData[,D]~as.numeric(MyData$H))
  
  #fit linear regression
  y=MyData[,D]
  x=MyData$H
  
  fit1 <- lm(y~poly(x,1,raw=TRUE))
  p1=summary(fit1)$coefficient[2,4]
  #if (p<=0.05){L=1}else{L=2}
  
  # fit the 2rd quadrilles regressions
  fit2 <- lm(y~poly(x,2,raw=TRUE))
  
  xx <- seq(0,100, length=50)
  summary(fit2)
  p2=summary(fit2)$coefficient[3,4]
  #if (p2<=0.05){L2=1}else{L2=2}
  if (p2<=0.001){p2="<0.001"}else{p2=p2}
  #lines(xx, predict(fit2, data.frame(x=xx)), col="blue",lty=L2)
  lines(xx, predict(fit2, data.frame(x=xx)), col="black",lty=1)
  
  if (p1<=0.001){p1="<0.001"}else{p1=p1}
  
  State = vector('expression',2)
  
  if (p1<=0.05){State[1] = substitute(expression(bold(p1:p-value == MYVALUE)), list(MYVALUE = format(p1, digits = 2)))[2]}else {State[1] = substitute(expression(p1:p-value == MYVALUE), list(MYVALUE = format(p1, digits = 2)))[2]}
  if (p2<=0.05){State[2] = substitute(expression(bold(p2:p-value == MYVALUE)),list(MYVALUE = format(p2, digits = 2)))[2]}else{State[2] = substitute(expression(p2:p-value == MYVALUE),list(MYVALUE = format(p2, digits = 2)))[2]}
  
  
  lines(xx, predict(fit1, data.frame(x=xx)), col="black",lty=1)
  
  legend("topright",legend=State, text.col=c("black","black"), bty = 'n',cex=0.8)
  
  axis(1, at = c(0,25,50,75,100), 
       label = rep("", 5), 
       tck = -0.04)
  axis(1, at = c(0,25,50,75,100), 
       label = c("0%","25%","50%","75%","100%"), 
       cex.axis=0.6,line = -0.7, lwd = 0)
  axis(1, at = 5, 
       label = "", 
       tck = -0.03)
  axis(1, at = 5, 
       label = "5%", 
       cex.axis=0.6,line = -1.1, lwd = 0)
}
###

Each_plot=function(DDay10,D,Xaxis,Yaxis,a,b,c){
  
  plot(as.numeric(DDay10$H),DDay10[,D], xaxt="n", yaxt="n",
       ylim=c(a,b),
       xlim=c(-5,105))
  
  # fit the 2rd quadrilles regressions
  y=DDay10[,D]
  x=DDay10$H
  fit2 <- lm(y~poly(x,2,raw=TRUE))
  xx <- seq(0,100, length=50)
  summary(fit2)
  p2=summary(fit2)$coefficient[3,4]
  if (p2<=0.05){L2=1}else{L2=2}
  if (p2<=0.001){p2="<0.001"}else{p2=p2}
  lines(xx, predict(fit2, data.frame(x=xx)), col="blue",lty=L2)
  
  # fit the linear regression
  l=lm(DDay10[,D]~as.numeric(DDay10$H))
  p=summary(l)$coefficient[2,4]
  if (p<=0.05){L=1}else{L=2}
  if (p<=0.001){p="<0.001"}else{p=p}
  abline(l,lty=L)
  State = vector('expression',2)
  State[1] = substitute(expression(italic(p1:p-value) == MYVALUE), 
                        list(MYVALUE = format(p, digits = 2)))[2]
  State[2] = substitute(expression(italic(p2:p-value) == MYVALUE), 
                        list(MYVALUE = format(p2, digits = 2)))[2]
  
  legend("topright",legend=State, bty = 'n',cex=0.8,text.col=c("black","blue"))
  #text(-4.5,b*0.99, label = Label1, cex=0.9)
  
  if (Xaxis==TRUE) {
    axis(1, at = c(0,25,50,75,100), 
         label = rep("", 5), 
         tck = -0.04)
    axis(1, at = c(0,25,50,75,100), 
         label = c("0%","25%","50%","75%","100%"), 
         cex.axis=0.6,line = -0.7, lwd = 0)
    axis(1, at = 5, 
         label = "", 
         tck = -0.03)
    axis(1, at = 5, 
         label = "5%", 
         cex.axis=0.6,line = -1.1, lwd = 0)
  }
  
  if (Yaxis==TRUE) {
    axis(2, tck = -0.04, at = seq(a,b, by = c), 
         label = rep("", length(seq(a,b, by = c))))
    axis(2, las=1, cex.axis=0.7, at = seq(a,b, by = c),
         line = -0.4, lwd = 0)
  }
  
}  

####
FigureS7_S8_Diveristy_Table=function(Table_old){
  #########Rarefaction######
  ps.rarefied = rarefy_even_depth(otu_table(t(Table_old),taxa_are_rows=F), rngseed=1, sample.size=2000, replace=T)
  
  Table_rarefied=t(ps.rarefied)
  Sample_all=colnames(Table_rarefied)
  Table=Table_rarefied[,which(sapply(strsplit(Sample_all, "H" , fixed = TRUE), "[", 1)=="CE1")]
  ##########################
  
  Richness=colSums(ifelse(Table> 1,1,0))
  Shannon=diversity(t(Table),"shannon")
  Simpson=diversity(t(Table),"simpson")
  invSimpson=diversity(t(Table),"invsimpson")
  Evenness=invSimpson/Richness
  Day=sapply(strsplit(sapply(strsplit(names(Richness), "D" , fixed = TRUE), "[", 2),"R",fixed=T),"[",1)
  H=sapply(strsplit(sapply(strsplit(names(Richness), "H" , fixed = TRUE), "[", 2),"D",fixed=T),"[",1)
  MyData=data.frame(Richness,Shannon,Simpson,Evenness,Day,H)
  return(MyData)
  ##########################
}

####
FigureS7_S8_10_12=function(MyData){
  Ordered_MyData=MyData[order(as.numeric(MyData$H)),]
  MyData=Ordered_MyData
  MyData$H=as.numeric(MyData$H)
  
  DDay10=filter(MyData, MyData$Day==10)
  DDay12=filter(MyData, MyData$Day==12)
  LIST=list(DDay10,DDay12)
  return(LIST)
}


