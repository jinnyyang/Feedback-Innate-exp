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
CorelationCoefficientEachH=function(Day,CE1_Only){
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
  plot(as.numeric(DDay12$H),DDay12[,D],ylim=c(a,b),xaxt="n", 
       xlim=c(-5,105))
  l=lm(DDay12[,D]~as.numeric(DDay12$H))
  
  #fit linear regression
  y=DDay12[,D]
  x=DDay12$H
  
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












DiveristyPlot_NoYaxis=function(x,y,a,b,c,Xaxis){
  plot(x,y, xaxt="n", yaxt="n",
       ylim=c(a,b),
       xlim=c(-5,105))
  
  # fit the 2rd quadrilles regressions
  fit2 <- lm(y~poly(x,2,raw=TRUE))
  xx <- seq(0,100, length=50)
  summary(fit2)
  p2=summary(fit2)$coefficient[3,4]
  if (p2<=0.05){L2=1}else{L2=2}
  if (p2<=0.001){p2="<0.001"}else{p2=p2}
  lines(xx, predict(fit2, data.frame(x=xx)), col="blue",lty=L2)
  
  # fit the linear regression
  l=lm(y~x)
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
}


DiveristyPlot_WithYaxis=function(x,y,a,b,c,Xaxis){
  plot(x,y, xaxt="n", yaxt="n",
       ylim=c(a,b),
       xlim=c(-5,105))
  axis(2, tck = -0.04, at = seq(a,b, by = c), 
       label = rep("", length(seq(a,b, by = c))))
  axis(2, las=1, cex.axis=0.7, at = seq(a,b, by = c),
       line = -0.4, lwd = 0)
  
  # fit the 2rd log regressions
  fit2 <- lm(y~poly(x,2,raw=TRUE))
  xx <- seq(0,100, length=50)
  Table_2=summary(fit2)
  p2=summary(fit2)$coefficient[3,4]
  if (p2<=0.05){L2=1}else{L2=2}
  if (p2<=0.001){p2="<0.001"}else{p2=p2}
  lines(xx, predict(fit2, data.frame(x=xx)), col="blue",lty=L2)
  
  # fit the linear regression
  l=lm(y~x)
  Table_1=summary(l)
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
  Table=list(Table_1,Table_2)
  return(Table)
}


##

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
