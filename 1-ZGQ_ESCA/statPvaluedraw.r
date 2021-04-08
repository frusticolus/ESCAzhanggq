# 针对用药前后进行细胞因子变化的显著性分析
# 先将数据变成短数据结构
library(tidyr)
library(readxl)
library(ggplot2)
library(ggpubr)

mydata <- read_excel("./2-analysis/zgq_cytokins.xlsx", sheet = 2)
widedt <- mydata
longdt <- gather(widedt, cytokins, measurement, `IL-1b`:VEGF, factor_key=TRUE)

# 所有的细胞因子名字
pdf("plots_pvalue.pdf")
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
    
    cydf2 <- subset(cydf,Path!="NA")
    p <- ggboxplot(cydf2, x = "Path", y = "measurement",
                   color = "Type", palette = "jco",
                   add = "jitter")
    pe <- p + ylab(cydf2$cytokins[1])+ stat_compare_means(aes(group = Type), method = "t.test")
    print(pe)
    
  }
}
dev.off()

# 分面
p <- ggboxplot(cydf2, x = "Type", y = "measurement",
               color = "Type", palette = "jco",
               add = "jitter",
               facet.by = "Path", short.panel.labs = FALSE)

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format")
