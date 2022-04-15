############### boxplot gene by gene

#########matrix type ######

###Single matrix, containing features/genes in columns and samples in rows, last column contain labels that is Normal, stage1, late stage ###
#########



df <-read.csv("3HCC_ext_val_mat", header=T, sep=",") 

require(reshape2)
df.m <- melt(df, id.var = "Gene")

require(ggplot2)
p <- ggplot(data = df.m, aes(x=variable, y=value)) 
p <- p + geom_boxplot(aes(fill = Gene)) + scale_fill_manual(values=c("pink", "cyan"))
# if you want color for points replace group with colour=Label
p <- p + geom_point(aes(y=value, group=Gene), position = position_dodge(width=0.75))
p <- p + facet_wrap( ~ variable, scales="free")
#p <- p + xlab("Genes") + ylab("Expression in RMA value") + ggtitle("The expression pattern of 17 signature genes in metastatic and primary tumor samples")
p <- p + xlab("Genes") + ylab("Expression Values)")
p <- p + guides(fill=guide_legend(title="Legend_Title"))
p 

