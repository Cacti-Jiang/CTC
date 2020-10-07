library(ggplot2)
library(gridExtra)
library(stringr) 

mat = read.table("CTC.exp.txt",head=T)
anno=read.table("anno_gene.txt",head=F,row.names=1)

mat <- mat[rownames(mat)%in%rownames(anno), ]
mat <- mat[rownames(anno), ]

G1 <- mat[rownames(anno)[anno$V2=="G1-S"], ]
G2 <- mat[rownames(anno)[anno$V2=="G2-M"], ]

G1_E <- log2(G1 + 1)
G2_E <- log2(G2 + 1)

G1_E_mean <- apply(G1[,], 1, function(x) log2(mean(x)+1))
G2_E_mean <- apply(G2[,], 1, function(x) log2(mean(x)+1))

G1_Er <- G1_E - G1_E_mean
G2_Er <- G2_E - G2_E_mean

E1 <- apply(G1_Er[,], 2, mean)
E2 <- apply(G2_Er[,], 2, mean)

score <- cbind(E1,E2)
E <- apply(score, 1, max)
score <- cbind(score, E)
score <- data.frame(score, pos=str_extract(colnames(mat), "HV|PA|PV|PoV"))

##info from heatmap
cycle <- read.table("cyc.cell.txt", header = F)

for(i in 1:nrow(score)){
    if(as.character(rownames(score)[i])%in%cycle$V1){
      score$heat[i] <- "Cycling"
    }else{
         score$heat[i] <- "Non-cycling"
      }
}


score$pos <- factor(score$pos, levels = c("HV", "PA", "PV", "PoV"))
p <- ggplot(score) + 
  geom_point(aes(x = E1, y = E2, colour = heat), size = 2) + 
  #  facet_grid(.~pos) 
  facet_grid(pos~.) 
p <- p + labs(x = "G1/S score", y = "G2/M score", title = "Cycling and Noncylcing CTC Cells Distribution 
              ",colour = "heat")+ theme_bw() +
  theme(panel.grid =element_blank(),plot.title = element_text(face = "bold", size=12, hjust = 0.5, lineheight=0.2))
p <- p + scale_shape_manual(values = 19) 
#  + scale_colour_brewer(palette = "Set1") 
p <- p + scale_colour_manual(name = "Cell Type", values=c("red",  "grey")) 

pdf(file = "pos.cyc.pdf", 6, 4)
p 
dev.off()
