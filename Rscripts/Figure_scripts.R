###Make a pie
library(ggplot2)
library(dplyr)
data <- data.frame(value=c(407,6157),
                   Gene= c("Accessory","Core" ))
#Compute percentages
data$fraction <- data$value / sum(data$value)
#Compute cumulative percentages (top of rectangle)
data$ymax <- cumsum(data$fraction)
#Compute bottom of each rectangle
data$ymin <- c(0, head(data$ymax, n=-1))
#Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2
#Compute a good label
data$label <- paste0(data$Gene, "\n ", data$value)
#make plot
g <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Gene))
p<- g +
  geom_rect()+
  geom_text(x=1.5, aes(y=labelPosition, label=label, color=Gene), size=5) + 
  #x here controls label position
  scale_fill_manual(values=c("#ce4257","#720026"))+
  scale_color_manual(values=c("#ce4257","#720026"))+
  coord_polar(theta="y")+
  xlim(c(-1,4))+
  theme_void()+
  theme(legend.position="none")
#export as 3x3in

###Histogram
fig2data<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/Fig2tbl.csv", sep=',', header=TRUE, skipNul = TRUE)
dat5 <- fig2data %>% select(c("Class", val="Count"))
g <- ggplot(dat5, aes(x=val, fill=Class))
b <- g + 
  geom_histogram(show.legend = FALSE,
                 binwidth=20,
                 boundary = 0)+
  labs(x = "Genomes",
       y = "Number of orthologs")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle=0, hjust=0.5, vjust=0),
        axis.text.y=element_text(size=12),
        axis.ticks = element_blank(),
        axis.line= element_blank(),
        axis.title = element_text(size = 15),
        panel.grid.minor=element_blank())+
  scale_y_continuous(limits=c(0,6500),breaks = seq(0, 8000, by = 2000))+
  scale_x_continuous(limits=c(-0,400),breaks = seq(0, 400, by = 100))+
  scale_fill_manual(values=c("#ce4257","#720026"))
#export as 6x4 in

###violin plot for orthogroups
extradata<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/Figure4tbl.csv", sep=',')
dat4 <- extradata %>% select(c("Strain", val="accessory_count"))
p <- ggplot(dat4, aes(x="Accessory", y=val))
v1<- p+
  geom_violin(fill="#ce4257", adjust=1)+
  theme_bw()+
  geom_boxplot(width=0.1)+
  scale_y_continuous(limits=c(175,275),breaks = seq(175,275, by = 25))+
  labs(x="",
       y="")+
  theme(panel.grid.minor=element_blank(),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=15))
p <- ggplot(extradata, aes(x="Core", y=core_count))
v2<- p+
  geom_violin(fill="#720026")+
  theme_bw()+
  geom_boxplot(width=0.1)+
  scale_y_continuous(limits=c(7100,7300),breaks = seq(6000, 7300, by = 50))+
  labs(x="",
       y="Count")+
  theme(panel.grid.minor=element_blank(),
        axis.title = element_text(size = 15),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=15))
library("ggpubr")
figure<- ggarrange(v2, v1, ncol=2)
#export as 4x6 in

###mean accessory genes with error bars
dat7<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/accessory_pergenome.csv", sep=',', header=TRUE)
p <- ggplot(dat7, aes(x=X, y=roll_ave, color="Accessory"))
l <- p + 
  geom_line()+
  geom_errorbar(aes(ymin=roll_ave-roll_st, ymax=roll_ave+roll_st), width=1, color="#ce4257", alpha=0.6)+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        panel.grid.minor=element_blank())+
  scale_color_manual(name="Gene Families",
                     values=c("#ce4257"),
                     guide="none")+
  scale_y_continuous(limits=c(150,450),breaks = seq(150, 450, by = 100))+
  scale_x_continuous(breaks = seq(0, 400, by = 100))+
  labs(x = "Genomes",
       y = "Accessory orthologs")
#export as 4x4.5 in

###Box plot for accuracy data
accuracy_dat<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/accuracy.csv", sep=',', header=TRUE) %>% select(c("Class","specificity","accuracy","recall"))
accuracy_dat <- accuracy_dat %>% select(c("Class", "specificity","recall","accuracy"))
accuracy_dat_melt<- melt(accuracy_dat, id.vars = "Class")
accuracy_dat_melt$Class <- factor(accuracy_dat$Class , levels=c("core","accessory"))
ggplot(accuracy_dat_melt,aes(y=value,x=variable,fill=Class))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        panel.grid.minor=element_blank())+
  scale_fill_manual(name="Ortholog class",values=c("#720026","#ce4257"))+
  labs(x = "Statistic",
       y = "Rate")
#export as 5x6in

###Phylogenetic tree with accessory genes
library(ggtree)
library(ape)
library(ggnewscale)
library(ggtreeExtra)
library(ggplot2)

tree <- read.tree("C:/Users/clair/OneDrive/Documents/Billmyre/DesjardGR17_Fig1_crypto_387_snpphylo.nwk")
droptips <- c("Bt92", "PLR10", "Bt206","WM626","H99")
tree_removedtips <- drop.tip(tree, droptips)

rooted_tree<-root(tree_removedtips,node=487,resolve.root = TRUE)
groupedtree<- groupClade(rooted_tree, .node=c(487,502,564))
groupedtree$edge.length[which(!(groupedtree$edge[,1] %in% groupedtree$edge[,2]))] <- sum(groupedtree$edge.length[which(!(groupedtree$edge[,1] %in% groupedtree$edge[,2]))])/2

extradata<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/Figure4tbl.csv", sep=',', header=TRUE)

dat1 <- extradata %>% select(c("Strain", "Group"))
dat2 <- extradata %>% select(c("Strain", "Mating.type","Clinical_vs_environmental"))
dat2<- melt(dat2,id="Strain")
dat4 <- extradata %>% select(c("Strain", val="accessory_count"))
#make tree with clade colors
p <- ggtree(tr=groupedtree, 
            layout= 'fan', 
            aes(color=group),
            open.angle=5,
            size=0.4)+
  geom_rootedge(rootedge=0.05, color="#B5838D",linewidth=0.4)
g<-p+
  geom_hilight(node=486,
               alpha=0.3, fill="grey",
               size=0.05)+
  geom_hilight(node=564,
               alpha=0.3, fill="grey",
               size=0.05)+
  geom_tiplab(color="grey",
              align=TRUE,
              linetype="dotted",
              linesize=0.1,
              show.legend=FALSE,
              size=1.3,
              offset=0.1,
              hjust=1,
              alpha=0.5)+
  geom_tiplab(color="black",
              align=TRUE,
              linetype=NULL,
              linesize=0.1,
              show.legend=FALSE,
              size=1.3,
              offset=0.1,
              hjust=1)+
  scale_colour_manual(
    name="Clade",
    values=c("#FFB4A2","#B5838D","#6D6875","#E5989B"),
    guide=guide_legend(order=1, override.aes = list(size=2, linewidth=2)),
    labels=c("VNII","VNI", "VNBII", "VNBI"))+
  new_scale_colour()
t<- g+
  geom_fruit(data=dat2, 
             geom=geom_tile,
             mapping=aes(y=Strain, x=variable, fill=value),
             color = "grey50", 
             offset = 0.18,
             pwidth=0.03)+
  scale_fill_manual(name="Metadata",
                    values=c("white","#332288","white","#44AA99","grey80"),
                    guide=guide_legend(order=2
                    ))+
  new_scale_fill()+
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7))
b<- t+
  geom_fruit(
    data=dat4,
    geom=geom_col,
    mapping=aes(y=Strain, x=val, fill="Accessory genes"),
    pwidth=0.4,
    axis.params = list(
      axis='x',
      text.size=1.5,
      nbreak=6,
      vjust=1,
      hjust=0.5,
      limits=c(0,300)),
    grid.params=list())+
  scale_fill_manual(
    name="Pangenome",
    values=c("#ce4257"))
#export as 10x12in

###make tree for long-read assemblies used
tree <- read.tree("C:/Users/clair/OneDrive/Documents/Billmyre/DesjardGR17_Fig1_crypto_387_snpphylo.nwk")
keeptips <- c("Muc525-1","PMHc1049.THER1.STOR","Gb159-1",
              "Bt105","PMHc1011.ENR","PMHc1023.ENR",
              "Bt35","PMHc1030.ENR.CLIN.ISO","Bt103","H99","Ftc555-1",
              "Bt63","Bt1","AD1_7a",
              "Bt22","Ftc239-1","T4",
              "125.91","Bt210","LP-RSA2296",
              "Bt130","A1-35-8",
              "A5-35-17","Bt106","Bt88",
              "C45","Ftc146-1",
              "Ze90_1","C23","C27","IUM-96-28-28","Bt133",
              "Bt65","Bt81","Bt89")
tree_removedtips <- keep.tip(tree, keeptips)
treeroot <- root(tree_removedtips, c("C45","T4","PMHc1023.ENR"),resolve.root = TRUE)
groupedtree<- groupClade(treeroot, .node=c(46,49,59))
groupedtree$edge.length[which(!(groupedtree$edge[,1] %in% groupedtree$edge[,2]))] <- sum(groupedtree$edge.length[which(!(groupedtree$edge[,1] %in% groupedtree$edge[,2]))])/2
g<-ggtree(groupedtree,layout="circular",aes(color=group),size=0.9)+
  geom_rootedge(rootedge=0.05, color="#B5838D",linewidth=0.9)+
  geom_tiplab(size=2.5,show.legend=FALSE,align=TRUE)+
  #geom_text(aes(label=node), size=2, vjust=0.5,color="black")+
  scale_colour_manual(
    name="Clade",
    values=c("#FFB4A2","#B5838D","#6D6875","#E5989B"),
    guide=guide_legend(order=1, override.aes = list(size=1, linewidth=2)),
    labels=c("VNI","VNII", "VNBII", "VNBI"))+
  scale_x_continuous(limits=c(-0.05,0.4))+
  theme(plot.margin = margin(t=0.5,r=0.5,l=1,b=1.5, "cm"))
#export as 5x7

###plot for accessory gene locations in H99
dat8<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/Orthos_H99.csv", sep=',', header=FALSE)
dat8<- subset(dat8, V2 != "MT")
H99_chr_length<-data.frame(chr_name=1:14,                        chr_len=c(2291499,1621675,1575141,1084805,1814975, 1422463,1399503,1398693,1186808,1059964,1561994,774062,756744,942867),
                           cen_start=c(978064,852348,1378287,708803,1564073,784407,525713,458008,803599,199433,870137,139632,581910,445442),
                           cen_end=c(1002059,892713,1409633,748322,1580221,819723,568089,505242,833856,230150,929453,170630,628588,471749))
mat_locus<-data.frame(chr_name=5,mat_start=158220,mat_end=329273)
ggplot(data=H99_chr_length)+
  labs(y="Base pair position",
       x="H99 chromosomes")+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12))+
  geom_segment(data=H99_chr_length,aes(x=chr_name,y=cen_start,yend=cen_end),
               color="grey90",
               size=5,
               lineend="round")+
  geom_segment(data=mat_locus,aes(x=chr_name,y=mat_start,yend=mat_end),
               color="lightblue",
               size=5,
               lineend="round")+
  geom_segment(data=dat8 %>%filter(V6=="accessory"), 
               aes(x=V2,xend=V2,y=V3,yend=V4),
               color="#ce4257",
               size=5)+
  geom_bar(aes(x=chr_name,y=chr_len),
           stat='identity',
           fill=NA,
           color='grey50',
           width=.7,
           linewidth=0.4)+
  coord_flip()+
  scale_x_reverse(breaks=seq(1,14, by=1))
#export as 4x6in

###make plot to show difference between accessory gene locations
library(ggsignif)
subtelo<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/longread_subtelo_acc.csv", sep=',', header=TRUE) %>% select(c("Genome","acc_per_bp_subtelo", "acc_per_bp_genome"))
subtelo_melt<- melt(subtelo, id.vars = "Genome")
subtelo_melt$variable <- factor(subtelo_melt$variable, levels=c("acc_per_bp_genome","acc_per_bp_subtelo"))
ggplot(subtelo_melt,aes(y=value*100000,x=variable))+
  geom_boxplot(fill="#ce4257")+
  geom_line(aes(group = Genome),color="grey40",alpha=0.3)+
  geom_point(position="identity")+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        panel.grid.minor=element_blank())+
  #scale_fill_manual(name="Ortholog class",values=c("#720026","#ce4257"))+
  labs(x = "Gene location",
       y = "Accessory genes per 100 kb")+
  scale_x_discrete(labels = c("Whole genome", "Subtelomeres"))+
  geom_signif(comparisons = list(c("acc_per_bp_genome", "acc_per_bp_subtelo")), 
              map_signif_level=TRUE)
#export as 4x4

###Chromosome 5 left arm genes
library(gggenomes)
s0<-tibble::tibble(
  seq_id=c("H99","125.91","T4","Bt1","Ftc146-1"),
  length=c(37000,37000,37000,37000,37000),
  start=c(2500,2500,2908,21654,31928),
  end=c(37000,37000,37000,37000,37000),)
g1 <- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/ident_island_genes.csv")
l2 <-read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/ident_island_links.csv")
p <- gggenomes(seqs=s0,
               genes =g1,
               links=l2,
               adjacent_only = FALSE) |> shift(all_of(s0$seq_id), by=s0$start)
p+
  geom_link(aes(fill=class),offset=0.15,color="NA")+
  geom_seq(size=0.5)+
  geom_seq_break(data_start = seqs(start < 2501),
                 data_end = seqs(end > 36000),
                 label="//")+
  geom_seq_label(size=3)+
  geom_gene(aes(fill=class))+
  geom_gene_tag(aes(label=feat_id),size=2,nudge_y=0.1)+
  scale_fill_manual(
    name="Gene class",
    values=c("#ce4257","#720026"))

###Create heatmap of ortholog conservation
library(pheatmap)
dat20g<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/one-to-one_orthocompare_groups.csv", sep=',', header=TRUE)
rownames(dat20g) <- dat20g[,1]
dat20g <- dat20g[,-1]
labels <- data.frame(Group = c("VNI", "VNII", "VNBI","VNBII"))
row.names(labels)<-c("VNI", "VNII", "VNBI","VNBII")
h<-pheatmap(dat20g,
            annotation_row = labels, 
            annotation_col = labels,
            show_rownames = FALSE,
            show_colnames = FALSE,
            color = hcl.colors(10, palette = "Rocket", alpha = TRUE, rev = TRUE, fixup = TRUE),
            silent = TRUE,
            annotation_colors=list(Group=c(VNI="#B5838D",VNII="#FFB4A2",VNBII="#6D6875",VNBI="#E5989B")))
#export as 6.5x5.5

###Create upset plot
library(ggupset)
mol_prop<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/percent_moltype.csv", sep=',', header=TRUE)%>% select(c("orthogroup","class","count","VNI_percent", "VNII_percent","VNBI_percent","VNBII_percent","mol_list"))
mol_specific<-mol_prop[mol_prop[,'VNI_percent']==0 | mol_prop[,'VNII_percent']==0 | mol_prop[,'VNBI_percent']==0 | mol_prop[,'VNBII_percent']==0 | mol_prop[,'class']=='accessory',]mol_specific$mol_list = strsplit(mol_specific$mol_list,",")
ggplot(mol_specific,aes(x=mol_list))+
  geom_bar(aes(fill=class))+
  scale_x_upset(name="Molecular type presence")+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, size=4)+
  scale_y_continuous(lim = c(0, 220), name = "Ortholog count")+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        panel.grid.minor=element_blank(),
        panel.grid.major.x=element_blank())+
  scale_fill_manual(
    name="Gene class",
    values=c("#ce4257","#720026"))
#save as 4x7.5

###Create plot showing accessory genes per genome curve for each molecular type
dat26<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/VNBII_accessory_pergenome.csv", sep=',', header=TRUE)
#dat27<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/VNII_accessory_pergenome.csv", sep=',', header=TRUE)
dat28<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/VNBI_accessory_pergenome.csv", sep=',', header=TRUE)
dat29<- read.csv("C:/Users/clair/OneDrive/Documents/Billmyre/VNI_accessory_pergenome.csv", sep=',', header=TRUE)
p <- ggplot()+ 
  geom_line(data=dat26,aes(x=X, y=roll_ave), color="#6D6875")+
  geom_errorbar(data=dat26,aes(x=X,ymin=roll_ave-roll_st, ymax=roll_ave+roll_st), width=1,
                color="#6D6875", alpha=0.6)+
  geom_line(data=dat28, aes(x=X, y=roll_ave), color="#E5989B")+
  geom_errorbar(data=dat28,aes(x=X,ymin=roll_ave-roll_st, ymax=roll_ave+roll_st), width=1, color="#E5989B", alpha=0.6)+
  geom_line(data=dat29, aes(x=X, y=roll_ave), color="#B5838D")+
  geom_errorbar(data=dat29,aes(x=X,ymin=roll_ave-roll_st, ymax=roll_ave+roll_st), width=1, color="#B5838D", alpha=0.6)+
  scale_y_continuous(limits=c(120,320),breaks = seq(120,320, by = 50))+
  scale_x_continuous(limits=c(20,200),breaks = seq(20, 200, by =60))+
  labs(x = "Genomes",
       y = "Accessory orthologs")+
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        panel.grid.minor=element_blank())
print(p)
#export as 4x5in

###Example of fitting power equations to data from accessory genes per genome
VNBI_rollave<- data.frame(x=dat28$X+1, y=dat28$roll_ave)
model <- model <- nls(y ~ a * x^b, data = VNBI_rollave, start = list(a = 1, b = 1))
eq <- substitute(italic(P) == a * italic(N)^b, 
                 list(a = format(unname(coef(model)[1]), digits = 5),
                      b = format(unname(coef(model)[2]), digits = 5)))
melt_df<-melt(dat28%>% select(-c(roll_ave, roll_st)),id="X")
curve_data <- data.frame(x = seq(0, max(VNBI_rollave$x)))
curve_data$y <- predict(model, newdata = curve_data)
p <- ggplot()+ 
  geom_point(data=melt_df,aes(x=X+1,y=value),alpha=0.3,color="grey",size=0.2)+
  geom_point(data=VNBI_rollave,aes(x=x,y=y), color="darkblue",size=1,shape=17)+
  geom_line(data = curve_data, aes(x = x, y = y), color = "red", linetype = "dashed",size=1.5) +
  theme_bw()+
  labs(x = "VNBI genomes",
       y = "Accessory orthologs")+
  theme(axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        legend.title=element_text(size=15),
        axis.text=element_text(size=12),
        panel.grid.minor=element_blank())+
  scale_x_continuous(limits=c(0,150),breaks = seq(0, 150, by =50))+
  annotate("text",x = 40, y = 10,label = as.expression(eq),size=4)
print(p)
#export as 4x5.5

###Accessory vs. species tree heatmap
tree <- read.tree("C:/Users/clair/OneDrive/Documents/Billmyre/Fig1_Crypto_Kwoniella_RelTime_timetree.nwk")
droptips <- c("'Cryp_sp._OR849'", "'Kwon_sp._B9012'","'Cryp_sp._OR918'", "'Kwon_sp._CBS9459'","'Kwon_sp._CBS6097'","'Kwon_sp._DSM27419'", "Trem_mese_ATCC28783","Cryp_gatt_MF34")
tree_removedtips <- drop.tip(tree, droptips)
rooted_tree<-root(tree_removedtips,outgroup=c("K_shandongensis_CBS12478","K_newhampshirensis_CBS13917"),resolve.root=TRUE)
rooted_tree$edge.length <- rooted_tree$edge.length * 0.25
rooted_tree$edge.length[which(!(rooted_tree$edge[,1] %in% rooted_tree$edge[,2]))] <- sum(rooted_tree$edge.length[which(!(rooted_tree$edge[,1] %in% rooted_tree$edge[,2]))])/2
species_ids<-read.csv("C:/Users/clair/Downloads/species_miniprot_identities.csv", sep=',', header=TRUE)
melt_species_ids <- melt(species_ids)
p <- ggtree(tr=rooted_tree,
            linewidth=0.2,
            layout = 'rectangular')+
  geom_rootedge(rootedge=2,linewidth=0.2)+
  geom_tiplab(size=2,offset=1)
r<-ggtree::rotate(p,26)
r2<-ggtree::rotate(r,27)
h<- r2+
  geom_fruit(data=melt_species_ids,
             aes(x=variable,y=X,fill=value),
             geom=geom_tile,
             offset=0.8,
             pwidth=6)+
  scale_fill_gradient(low = "white", high = "blue",na.value="white",
                      name="Identity")+
  theme(legend.key.size = unit(0.1, "in"))
print(h)
#export as 4.5x 11in

###Gene tree
library(gggenomes)
s0<-tibble::tibble(
  bin_id=c("Chr14", "ChrH","Chr04","ChrI"),
  seq_id=c("Chr14","ChrH", "Chr04","ChrI"),
  length=c(967068,1265488,1053706,989306),
  start=c(1,1,867436,830932),
  end=c(126713,126152,1053706,989306)
)
l1 <-read.csv("C:/Users/clair/Downloads/CN06540_links.csv")
p<-gggenomes(seqs=s0,links=l1)
p2<-
  p + 
  geom_seq_label()+
  geom_link_curved(fill="blue", color="blue",alpha=0.5, color="NA",offset=0.0,
                   curve=6)+
  geom_seq(size=1,fill="NA") +
  geom_seq_break(size=6)
p3<- p2 %>% flip(3,4)