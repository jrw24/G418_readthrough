library("ggplot2")
library("plyr")
library("reshape2")
library("scales")

args= commandArgs(trailingOnly = TRUE)
fname= args[1]
wd= args[2]
setwd(args[2])
output= args[3]

biglist= list()

pos= c(-50:149)
biglist[["pos"]]= pos

ftsize= c(15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35)
#ftsize= c(29,30,31,32,33)
#ftsize= c(18,19,20,21,22,23)

for (ft in ftsize) 
{
  infile= sprintf("%s_%sf_avg_1.csv",fname,ft)
  input= read.csv(infile, header= TRUE, sep= ",")
  attach(input)
  name= paste("a", ft, sep= "")
  biglist[[name]]= c(avg)
  detach(input)
}

dt= data.frame(biglist)
dtT= t(dt)

# Filter size of interest
dtF= dtT#[,30:80]
write.table(dtF,"biglist_F_start.csv", row.names= TRUE, col.names= FALSE, sep=",")
dtW= read.csv("biglist_F_start.csv", header= TRUE, sep= ",")

dt.m= melt(dtW)
dt2.m= ddply(dt.m, .(pos), transform, rescale= rescale(value))
#dt2.m= ddply(dt.m, .(pos), transform, value)

plot= ggplot(dt2.m, aes(variable, pos)) + 
  geom_tile(aes(fill= value), color= "white") + 
  scale_fill_gradient(low= "white", high= "black")+
  scale_x_discrete(name= "5'ends of reads to START", labels= c(-50:149))+
  scale_y_discrete(name= "Footprint size", labels= ftsize)

heatmap= plot+
  theme(panel.background= element_blank(),
        panel.grid= element_blank(), 
        panel.border= element_blank(),
        axis.text.x= element_text(size= 4, color='black'),
        axis.text.y= element_text(size= 6, color='black'),
        axis.ticks.x= element_line(size= 0.25, color='transparent'),
        axis.ticks.y= element_line(size= 0.25, color='transparent'),
        axis.title= element_text(size= 8, face= 'plain'),
        plot.margin= unit(c(3,3,3,3),"mm"),
        legend.position= "right",
        legend.title= element_blank(),
        legend.box.margin= margin(t= 0, r= 0, b= 0, l= 0, unit= "cm"),
        legend.key.size= unit(0.2, units= "cm"),
        legend.text= element_text(size= 3))

ggsave(file= output, heatmap, width= 6, height= 2, units='in', dpi= 600, pointsize= 0.5)
plot(heatmap)
dev.off()
