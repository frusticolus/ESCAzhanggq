library(data.table)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(grid)
library(cowplot)

setwd(getwd())
ab_result <- fread("data/CIBERSORTx_ESPD_absolute_Results.txt")
#ab_result <- fread("data/CIBERSORTx_ESPD_relative_Results.txt")


sort_cell <- sort(colnames(ab_result)[2:23])
class <- ifelse(sort_cell %like% "B cells", "B cell", sort_cell)
class <- ifelse(class %like% "T cells", "T cells", class)
class <- ifelse(class %like% "NK cells", "NK cells", class)
class <- ifelse(class %like% "Dendritic", "Dendritic", class)
class <- ifelse(class %like% "Mast", "Mast", class)
class <- ifelse(class %like% "Macrophages", "Macrophages", class)

raw_pal <- pal_npg("nrc")(10)[c(10, 1:3,5:9,4)]

create_pal <- function(class, raw_pal){
  tb <- table(class)
  if (length(raw_pal) < length(tb)){
    stop(paste0("length of raw pal should greater or eqal than ", length(tb)))
  }
  pal <- c()
  for (i in 1:length(tb)){
    tmp_pal <- raw_pal[i]
    pal <-  c(pal, alpha(tmp_pal, (seq(.4, 1, length.out = tb[i]))))
  }
  return(pal)
}

# give pal
pal <- create_pal(class, raw_pal)
names(pal) <- sort_cell

# transorm data
ab_dp <- melt(ab_result[, 1:23], id.vars = "Patient", measure.vars = colnames(ab_result)[2:23],
              variable.name = "CellType", value.name = "Cell Proportion")



# fix cell type order
ab_dp$CellType <- factor(ab_dp$CellType, levels = sort_cell)
# fix patient order 
# not run in relative mode !!!
order_patient <- ab_dp[,.(total_fraction = sum(`Cell Proportion`)),by=Patient][order(-total_fraction), Patient]

ab_dp <- ab_dp[order_patient,on="Patient"]
ab_dp$Patient <- factor(ab_dp$Patient, levels = unique(ab_dp$Patient))

clin = read.table("data/sample_clinical_info.xls", header = T, sep = "\t",stringsAsFactors = F)
head(clin)
clin = clin[which(clin$基线_RNA != 0),]
clin$Response_imaging=ifelse(clin$Response_imaging == "","NA",clin$Response_imaging)
clin$Response_imaging=factor(clin$Response_imaging,levels = c("PR","SD","PD"))
clin$Response_irPRC=ifelse(clin$Response_irPRC == "","NA",clin$Response_irPRC)
clin$Response_irPRC = factor(clin$Response_irPRC,levels = c("pCR","MPR","noMPR","NA"))
clin$Response1=as.character(clin$Response1)
clin$Response2=as.character(clin$Response2)
clin$Alias = factor(clin$Alias, levels = unique(ab_dp$Patient))

# plot figure
p1<-ggplot(ab_dp, aes(x = Patient, y = `Cell Proportion`, fill = CellType)) + geom_bar(stat = "identity", width = 0.8) + theme_bw() +
  scale_fill_manual(values = pal)  + labs(fill = "Cell type")+
  #facet_wrap(~ITH_level, nrow = 1, scales = "free_x")+ 
  labs(x="")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_blank()) +
  theme( axis.text.y=element_text(colour = "black"),
         axis.ticks.y=element_line(colour = "black")) + 
  theme( #axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.line.x = element_blank()) + 
  scale_y_continuous(expand = c(0, 0))


p2<-ggplot(data=clin,aes(x=Alias,y=1,fill=Response_irPRC))+
  geom_bar(stat = "identity", width=.9)+
  theme_void()+
  #guides(fill = guide_legend(title = 'irPRC'))+
  theme(panel.spacing.x = unit(1, "mm"))+
  scale_fill_manual("irPRC", values = c("pCR" = "#2e397d", "MPR" = "#00837c", "noMPR" = "#00cc3d","NA" = "lightgray"))

p3<-ggplot(data=clin,aes(x=Alias,y=1,fill=Response_imaging))+
  geom_bar(stat = "identity", width=.9)+
  theme_void()+
  guides(fill = guide_legend(title = 'Response_imaging'))+
  theme(panel.spacing.x = unit(1, "mm"))+
  scale_fill_manual("Response_imaging", values = c("PR" = "#2e397d", "SD" = "#00837c", "PD" = "#00cc3d"))

p4<-ggplot(data=clin,aes(x=Alias,y=1,fill=Response1))+
  geom_bar(stat = "identity", width=.9)+
  theme_void()+
  #guides(fill = guide_legend(title = 'Response1'))+
  theme(panel.spacing.x = unit(1, "mm"))+
  scale_fill_manual("Response1", values = c("1" = "#2e397d", "2" = "#00837c", "3" = "#00cc3d"))

p5<-ggplot(data=clin,aes(x=Alias,y=1,fill=Response2))+
  geom_bar(stat = "identity", width=.9)+
  theme_void()+
  #guides(fill = guide_legend(title = 'Response2'))+
  theme(panel.spacing.x = unit(1, "mm"))+
  scale_fill_manual("Response2", values = c("1" = "#2e397d", "2" = "#00837c", "3" = "#00cc3d","4"="#ffeb00"))

legend1 <- plot_grid( get_legend(p1), ncol = 1)
legend2 <- plot_grid( get_legend(p2), ncol = 1)
legend3 <- plot_grid( get_legend(p3), ncol = 1)
legend4 <- plot_grid( get_legend(p4), ncol = 1)
legend5 <- plot_grid( get_legend(p5), ncol = 1)
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")
p5 <- p5 + theme(legend.position = "none")

plot <- plot_grid(p1,p2,p3,p4,p5, align = "v", ncol = 1, axis = "tb",rel_heights  = c(3, 0.1,0.1,0.1,0.1))
legend_group <- plot_grid(legend2,legend3,legend4, legend5, hjust = 0, vjust = 1,scale = c(1., 1., 0.9, 0.9))
legend <- plot_grid(legend1,legend_group,ncol = 1,rel_heights=c(0.6, 0.4))
plot_grid(plot, legend, nrow = 1,rel_widths=c(0.6, 0.4))

ggsave("output/ESPD.abolute_bar.pdf", width = 9, height = 6)
#ggsave("output/ESPD.relative_bar.pdf", width = 9, height = 6)




