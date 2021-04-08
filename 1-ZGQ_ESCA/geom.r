library('readxl')
library('ggplot2')

mydata <- read_excel("./2-analysis/zgq_cytokins.xlsx", sheet = 2)

#header <- gsub("\\(\\d+\\)","",colnames(mydata))
#colnames(mydata) <- header

for (i in c(7:33)) {
  print(i)
  aa = length(which(mydata[,i]=='OOR<'))
  print(aa)
}

# 以细胞因子为维度，绘制单一细胞因子用药前后的变化
col <- colnames(mydata)
pdf("cytokins_plots.pdf")
for (i in 8:34) {
  ytxt <- col[i]
  mydf <- mydata[,c('Name','Type','Path', col[i])]
  aa <- as.numeric(unlist(mydf[,4]))
  p <- ggplot(data = mydf, mapping = aes(x=Type, y=aa, colour=Path, group=Name)) + 
    geom_line() + 
    ylab(ytxt) 
  print(p)
  #ggsave(p,file=paste(ytxt,".pdf",sep=""))
}
dev.off()

# 以患者为维度，绘制每个患者用药前后细胞因子变化差异图
pts <- unique(mydata$ID)
pdf("plots.pdf")
for (i in pts) {
  ## 选取所需的行
  mydf <- dplyr::filter(mydata, grepl(i, ID))
  
  ## 把含有OOR变成NA
  mydf[mydf=="OOR<"]<-NA

  ## 去除掉含有OOR的列
  mydf2 <- t(na.omit(t(mydf)))
  path <- mydf$Path[1]
  
  mydf2 <- t(mydf2)
  before <- as.numeric(mydf2[8:dim(mydf2)[1],][,1])
  drug <- as.numeric(mydf2[8:dim(mydf2)[1],][,2])
  delta <- drug - before
  path <- rep(path,length(delta))
  name <- rownames(mydf2)[8:dim(mydf2)[1]]
  out <- data.frame(cytokins=name, value=delta, response=path)
  
  if ('pCR' %in% mydf$Path) {
    col = 'steelblue'
  } else if ('MPR' %in% mydf$Path) {
    col = 'green'
  } else if ('noMPR' %in% mydf$Path) {
    col = 'red'
  } else {
    col = 'grey'
  }
  p<- ggplot(data = out, aes(x=value, y=cytokins)) + 
    geom_bar(stat="identity",fill=col) +
    theme_classic() + 
    labs(title=i)
  
  print(p)
}
dev.off()

# 以患者为维度，绘制用药前后的细胞因子数值，同时进行显著性分析
# 先将数据变成短数据结构

widedt <- mydata
longdt <- gather(widedt, cytokins, measurement, `IL-1b`:VEGF, factor_key=TRUE)

## 挑选细胞因子
ptdf <- dplyr::filter(longdt, grepl('Y4707497', ID))
## 把测量值变成数值型变量
ptdf$measurement <- as.numeric(ptdf$measurement)
## 如果是NA数值变成0
ptdf$measurement[is.na(ptdf$measurement)] <- 0

## 下面是不同的处理方法，针对OOR<的处理
#ptdf[ptdf=="OOR<"]<-'0'
#ptdf2 <- na.omit(ptdf)
#ptdf$measurement <- as.numeric(ptdf$measurement)

ggplot(ptdf, mapping = aes(x=cytokins, y=measurement, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


# 针对用药前后进行细胞因子变化的显著性分析
# 先将数据变成短数据结构
library(tidyr)
library('readxl')
library('ggplot2')
mydata <- read_excel("./2-analysis/zgq_cytokins.xlsx", sheet = 2)
widedt <- mydata
longdt <- gather(widedt, cytokins, measurement, `IL-1b`:VEGF, factor_key=TRUE)

# 所有的细胞因子名字
cynames <- colnames(mydata)[8:34]
for (i in cynames) {
  ## 挑选细胞因子
  cydf <- dplyr::filter(longdt, grepl(i, cytokins))
  ## 把测量值变成数值型变量
  cydf$measurement <- as.numeric(cydf$measurement)
  # 会出现NA行，需要把出现NA行的患者找到
  rmpts <- cydf[complete.cases(cydf)==FALSE, ]$ID
  cydf <- dplyr::filter(cydf,  !cydf$ID %in% rmpts)

  ## 计算用药前后是否存在显著差异
  pre <- dplyr::filter(cydf, grepl('Before', Type))
  che <- dplyr::filter(cydf, grepl('Drug', Type))
  prev <- pre$measurement
  chev <- che$measurement
  num <- length(prev)
  if (num < 2 ) {
    cat(i, "NA", num, sep="\t", "\n")
  } else {
    pvalue <- (t.test(prev,chev, paired=T))$p.value
    cat(i, pvalue, num, sep="\t", "\n")
  }
}


## 绘制用药前后某细胞因子的baxplot图
ggplot(cydf, mapping = aes(x=Type, y=measurement)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
