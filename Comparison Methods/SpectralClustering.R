install.packages('kernlab')
library(kernlab)
library(NMI);

setwd("D:/R")

sets<-c('2 whiteside','4 banknote-authentication','4 iris','5 userknowledge','10 BreastTissue','28 foresttypes training','71 SonyAIBORobotSurface1','31 wdbc','78 mice protein','81 ProximalPhalanxOutlineAgeGroup','129 SwedishLeaf','145 Plane');
k<-c(2,2,3,4,6,4,2,2,8,3,15,7);


sets<-c("2 whiteside","4 banknote-authentication","4 iris","10 BreastTissue","78 mice protein","129 SwedishLeaf");


for(i in 1:length(sets)){
  name = sets[i]
  name = paste(name,'Dip',sep='');
  
  #name = paste(name,'Dip',sep='')
  dat = paste(name,'.txt',sep='')
  data<-read.csv2(dat,header = FALSE,sep=",", dec=".")
  dim = dim(data)[2]
  
  mydata<-data[,-dim]
  da<-data[,dim]
  mydata <- as.matrix(mydata);
  
  kcluster=length(unique(da))
  print(name)
  print(kcluster)
  
  for(j in 1:100){
    sc <- specc(mydata, centers=kcluster)
    print(j)
    la<-sc@.Data
    typeof(la)
    
    namel = paste(name,'_spec',sep='')
    namel = paste(namel,j,sep='')
    namel = paste(namel,'.csv',sep='')
    namel
    write.table(la, file = namel,col.names=FALSE, row.names = FALSE,sep=", ");
  }
}
