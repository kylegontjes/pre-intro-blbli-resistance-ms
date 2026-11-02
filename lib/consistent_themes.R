# Consistent themes

## Packages
library(gridExtra)
library(ggtext)
library(kableExtra)
library(knitr)
library(hues)

# Global font: Arial (Nature: Helvetica/Arial)
theme_set(theme_get() + theme(text       = element_text(family = "Arial", colour = "black"),
                              axis.text  = element_text(colour = "black"),
                              axis.ticks = element_line(colour = "black")))
update_geom_defaults("text",  list(family = "Arial"))
update_geom_defaults("label", list(family = "Arial"))

# Manuscript Background theme
theme_bw_me <- theme(panel.background = element_rect(fill = "white",colour = NA), panel.grid = element_blank(),
                     strip.background = element_rect(fill = "white",colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.line = element_line(colour = "black"),legend.position = "bottom")

# Resistance categories
resistance_cat_colors <- c("Susceptible" = "#005AB5","Intermediate"="#FFC20A","Resistant" = "#DC3220")
resistance_cat_scale <- scale_fill_manual(values=resistance_cat_colors, name="Resistance category", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))
resistance_cat_scale_v <- scale_fill_manual(values=resistance_cat_colors, name="Resistance category", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))

# Resistance categories
feature_colors <- c(`1` = "black",`0`="white")
feature_scale <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Tip state", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))
feature_scale_v <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Tip state", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))
plasmid_scale_v <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Plasmid", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))
plasmid_scale_v_key <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Plasmid", guide = guide_legend(ncol=1, title.position = "top", label.position = "right",override.aes = list(colour = "black")))
genotype_scale_v <- scale_fill_manual(breaks = c(1,0), values = c('black','white') ,labels=c("Present","Absent"),name="Genotype", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))

MVB_IR_scale <- scale_fill_manual(breaks = c("IR_num_log_2_diff","MVB_num_log_2_diff"),values=c(5,6),labels = c("Imipenem-relebactam","Meropenem-vaborbactam"),name="Antibiotic", guide = guide_legend(ncol=1,title.position = "top", label.position = "right")) 

# Non-Susceptibility
resistance_prop_scale <- scale_fill_manual(breaks = c("blbli_res_prop", "blbli_sus_prop"),
                                           values = c("black", "white"),
                                           labels = c("Resistant", "Susceptible"),
                                           name = "Resistance",
                                           guide = guide_legend(ncol = 1, title.position = "top", label.position = "right",
                                                                override.aes = list(colour = "black", linewidth = 0.2)))
resistance_scale <- scale_fill_manual(breaks = c("Non-Susceptible","Susceptible"),values=c("black","white"),labels = c("Resistant","Susceptible"),name="Resistance",guide = guide_legend(order=3,title.position = "top", label.position = "right",nrow=2),drop = FALSE )
resistance_scale_key <- scale_fill_manual(breaks = c("Non-Susceptible","Susceptible"),values=c("black","white"),labels = c("Resistant","Susceptible"),name="Resistance",guide = guide_legend(order=3,title.position = "top", label.position = "right",nrow=2,override.aes = list(colour = "black")),drop = FALSE )
resistance_scale_order_2 <- scale_fill_manual(breaks = c("Non-Susceptible","Susceptible"),values=c("black","white"),labels = c("Resistant","Susceptible"),name="Resistance",guide = guide_legend(order=2,title.position = "top", label.position = "right",nrow=2),drop = FALSE )
resistance_scale_order_1 <- scale_fill_manual(breaks = c("Non-Susceptible","Susceptible"),values=c("black","white"),labels = c("Resistant","Susceptible"),name="Resistance",guide = guide_legend(order=1,title.position = "top", label.position = "right",nrow=2),drop = FALSE )

# MIC Numeric
# MIC Numeric
colfunc <- colorRampPalette(c("white", "red")) 
log_breaks <- c(-2,-1,0,1,2,3,4,5)
Log2_scale <-scale_fill_manual(breaks = log_breaks,values = colfunc(8),labels = c("≤0.25/IN"," 0.5/IN"," 1/IN"," 2/IN"," 4/IN"," 8/IN"," 16/IN","≥32/IN"),name = "Minimum inhibitory concentration (MIC)",guide = guide_legend(title.position = "top",label.position = "bottom",nrow=1,keywidth=2,order=4,override.aes = list(colour = "black", linewidth = 0.2)),drop = FALSE ,limits=force)   # CHANGED
Log2_scale_key <-scale_fill_manual(breaks = log_breaks,values = colfunc(8),labels = c("≤0.25/IN"," 0.5/IN"," 1/IN"," 2/IN"," 4/IN"," 8/IN"," 16/IN","≥32/IN"),name = "Minimum inhibitory concentration (MIC)",guide = guide_legend(title.position = "top",label.position = "bottom",nrow=1,keywidth=2,order=4,override.aes = list(colour = "black", linewidth = 0.2)),drop = FALSE ,limits=force)   # CHANGED: linewidth 0.2
Log2_scale_3_order <-scale_fill_manual(breaks = log_breaks,values = colfunc(8),labels = c("≤0.25/IN"," 0.5/IN"," 1/IN"," 2/IN"," 4/IN"," 8/IN"," 16/IN","≥32/IN"),name = "Minimum inhibitory concentration (MIC)",guide = guide_legend(title.position = "top",label.position = "bottom",nrow=1,keywidth=2,order=3,override.aes = list(colour = "black", linewidth = 0.2)),drop = FALSE ,limits=force)   # CHANGED
Log2_scale_gwas_figure <- scale_fill_manual(breaks = log_breaks, values = colfunc(8),
                                            labels = c("≤0.25/IN", " 0.5/IN", " 1/IN", " 2/IN", " 4/IN", " 8/IN", " 16/IN", "≥32/IN"),
                                            name = "Minimum inhibitory concentration (MIC)",
                                            guide = guide_legend(title.position = "top", label.position = "bottom", nrow = 1,
                                                                 keywidth = unit(5, "mm"), keyheight = unit(2.5, "mm"),
                                                                 override.aes = list(colour = "black", linewidth = 0.2),          # CHANGED
                                                                 theme = theme(legend.text = element_text(size = 6, hjust = 0.5,
                                                                                                          margin = margin(t = 2, unit = "pt"))),
                                                                 order = 4),
                                            drop = FALSE, limits = force)

Log2_scale_gwas_s_figure_3 <- scale_fill_manual(breaks = log_breaks, values = colfunc(8),
                                                labels = c("≤0.25/IN", "0.5/IN", "1/IN", "2/IN", "4/IN", "8/IN", "16/IN", "≥32/IN"),
                                                name = "Minimum inhibitory concentration (MIC)",
                                                guide = guide_legend(title.position = "top", label.position = "bottom", nrow = 1,
                                                                     keywidth = unit(5, "mm"), keyheight = unit(2.5, "mm"),
                                                                     override.aes = list(colour = "black", linewidth = 0.2),          # black key outline
                                                                     theme = theme(legend.text = element_markdown(size = 5.5, hjust = 0.5,
                                                                                                                  margin = margin(t = 2, unit = "pt"))),
                                                                     order = 4),
                                                drop = FALSE, limits = force)
Log2_scale_s6 <- scale_fill_manual(breaks = log_breaks, values = colfunc(8),
                                   labels = c("≤0.25/IN", "0.5/IN", "1/IN", "2/IN", "4/IN", "8/IN", "16/IN", "≥32/IN"),
                                   name = "Minimum inhibitory concentration (MIC)",
                                   guide = guide_legend(title.position = "top", label.position = "bottom", nrow = 1, order = 3,
                                                        keywidth = unit(3.5, "mm"), keyheight = unit(2.5, "mm"),
                                                        override.aes = list(colour = "black", linewidth = 0.2),
                                                        theme = theme(legend.text = element_text(size = 6, hjust = 0.5,
                                                                                                 margin = margin(t = 1, unit = "pt")),
                                                                      legend.key.spacing.x = unit(0.5, "mm"))),
                                   drop = FALSE, limits = force)

## X-axis
MIC_x_scale <- scale_x_continuous(breaks=log_breaks, labels=c("≤0.25","0.5", "1","2","4","8","16","≥32"))
## MIC Scales
MIC_IR_x_scale <- scale_x_continuous(breaks=c(-2,-1,0,1,2,3,4,5),
                   labels=c("≤0.25/4","0.5/4", "1/4","2/4","4/4","8/4","16/4","≥32/4"))
MIC_MVB_x_scale <- scale_x_continuous(breaks=c(-2,-1,0,1,2,3,4,5),
                                      labels=c("≤0.25/8","0.5/8", "1/8","2/8","4/8","8/8","16/8","≥32/8")) 

# Clade 
clade_colors <- c("Clade I" = "red","Clade II"="blue")
clade_colors_scale <- scale_fill_manual(values=clade_colors, name="Clade", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))
clade_colors_scale_v <- scale_fill_manual(values=clade_colors, name="Clade", guide = guide_legend(order=1,ncol=1, title.position = "top", label.position = "right"))
clade_colors_scale_3_col<- scale_fill_manual(values=clade_colors, name="Clade", guide = guide_legend(order=1,nrow=3, title.position = "top", label.position = "right"))
clade_colors_scale_point <- scale_color_manual(values=clade_colors, name="Clade", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))

### Clustering
cluster_colors <- c("No feature" = "white","Singleton" = "black","Cluster 1" = "green","Cluster 2" = "orange","Cluster 3" = "purple","Cluster 4" = "blue","Cluster 5" = "red","Cluster 6" = 'pink',"Cluster 7" = "#F0E442","Cluster 8" = '#00CED1')
cluster_labels <- names(cluster_colors) 
cluster_colors_gray <- c("No feature" = "gray","Singleton" = "black","Cluster 1" = "green","Cluster 2" = "orange","Cluster 3" = "purple","Cluster 4" = "blue","Cluster 5" = "red","Cluster 6" = 'pink',"Cluster 7" = "#F0E442","Cluster 8" = '#00CED1')
cluster_labels_gray <- names(cluster_colors) 
cluster_scale <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=3, title.position = "top", label.position = "right"))
cluster_scale_1_order <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=1,ncol=3, title.position = "top", label.position = "right"))
cluster_scale_1_order_2_col <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=1,ncol=2, title.position = "top", label.position = "right"))
cluster_scale_4 <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels  %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,nrow=4, title.position = "top", label.position = "right"))
cluster_scale_2_col <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels  %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=2, title.position = "top", label.position = "right"))
cluster_scale_2_col_key <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels  %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=2, title.position = "top", label.position = "right",override.aes = list(colour = "black")))
cluster_scale_3_col_key <- scale_fill_manual(values=cluster_colors,breaks = cluster_labels,labels = cluster_labels  %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=3, title.position = "top", label.position = "right",override.aes = list(colour = "black")))
cluster_scale_gray <- scale_fill_manual(values=cluster_colors_gray,breaks = cluster_labels_gray,labels = cluster_labels_gray %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=3, title.position = "top", label.position = "right"))
cluster_scale_gray_color <- scale_color_manual(values=cluster_colors_gray,breaks = cluster_labels_gray,labels = cluster_labels_gray %>% recode(.,"No feature" = "Susceptible"),name="Phylogenetics of resistance", guide = guide_legend(order=2,ncol=3, title.position = "top", label.position = "right"))

# Panels
panel_scale <- scale_fill_manual(breaks = c("known","novel_GWAS","known_ST258_GWAS"),values=hues::iwanthue(5),labels = c("Non-carbapenemase genotypes","Novel GWAS hits","Non-carbapenemase + novel GWAS"),name="Genotype panel",guide = guide_legend(ncol=3,title.position = "top", label.position = "right",order = 1)) 
panel_scale_v <- scale_fill_manual(breaks = c("known","novel_GWAS","known_ST258_GWAS"),values=hues::iwanthue(5),labels = c("Non-carbapenemase genotypes","Novel GWAS hits","Non-carbapenemase + novel GWAS"),name="Genotype panel",guide = guide_legend(ncol=1,title.position = "top", label.position = "right",order = 1)) 

# OmpK36 scale
ompK36_colors <- c("No loop 3 insertion or PFAV" = "black","Putative function-altering variant (PFAV)"="#00AFBB","Loop 3 insertion"="#FC4E07")
ompK36_scale <- scale_fill_manual(breaks = names(ompK36_colors),values=ompK36_colors,name="Status of the OmpK36 porin", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))
ompK36_color_scale <- scale_color_manual(breaks = names(ompK36_colors),values=ompK36_colors,name = "Status of the OmpK36 porin", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))

# Tn4401 allele
tn4401_alleles <- c("Tn4401a","Tn4401b","Tn4401d","Tn4401 del 1-3391 6920-7126","Tn4401 del 1-554 7008-7075","Tn4401 del 6920-7126")
tn4401_alleles_md <- gsub("Tn4401","<i>tn4401</i>",tn4401_alleles)
tn4401_scale <- scale_fill_manual(breaks=tn4401_alleles,labels = tn4401_alleles_md,values = c(8,9,10,11,12,13),name = "<i>tn4401</i> isoform",guide = guide_legend(order=5,title.position = "top", label.position = "right",ncol=1))

#KPC Scale
KPC_values <- c("blaKPC-2","blaKPC-3","blaKPC-5","blaKPC-10","No blaKPC")
KPC_gene_values <- gsub("bla","",KPC_values)
format_blaKPC_md <- function(gene) {
  if(grepl("^blaKPC-", gene)){
    allele <- sub("blaKPC-", "", gene)
    paste0("<i>bla</i><sub>KPC-", allele, "</sub>")
  }  else {
    if(grepl("No blaKPC",gene) == T){
      gsub("blaKPC","<i>bla</i><sub>KPC</sub>",gene)
    } else {
      gene
    }
  } 
}

KPC_names <- c(paste0('bla',KPC_gene_values %>% subset(.!='No KPC')),"No blaKPC")

KPC_values_md <- sapply(KPC_values,format_blaKPC_md) 

kpc_gene_scale <- scale_fill_manual(breaks=KPC_gene_values,labels = KPC_names,values = c(5,6,7,4,"white"),name = "blaKPC allele",guide = guide_legend(order=4,title.position = "top", label.position = "right",ncol=1))
 
kpc_scale <- scale_fill_manual(breaks=KPC_values,labels = KPC_values_md,values = c(5,6,7,4,"white"),name = "<i>bla</i><sub>KPC</sub> allele",guide = guide_legend(order=6,title.position = "top", label.position = "right",ncol=1))
kpc_scale_h <- scale_fill_manual(breaks=KPC_values,labels = KPC_values_md,values = c(5,6,7,4,"white"),name = "<i>bla</i><sub>KPC</sub> allele",guide = guide_legend(order=6,title.position = "top", label.position = "right",nrow=1))
 
# AA552
AA552_fill <-  scale_fill_manual(breaks = c("Present","Absent"),values=c("brown","gray"),name="AA552 <i>bla</i><sub>KPC</sub> plasmid", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))
AA552_fill_1_col <-  scale_fill_manual(breaks = c("Present","Absent"),values=c("brown","gray"),name="AA552 <i>bla</i><sub>KPC</sub> plasmid", guide = guide_legend(ncol=1, title.position = "top", label.position = "right"))


# Variant type
variant_type_scale <-  scale_fill_manual(breaks = c("SNP","INDEL","Insertion"),values = c("orange","purple","darkgray"), labels = c("SNP","INDEL","Insertion sequence"),name="Variant type", guide = guide_legend(nrow=1, title.position = "top", label.position = "right"))
 
## Mutation scale
mutation_scale <- scale_fill_manual(values = c("#b3a2c7",'#5e3c99','darkgray',"white"),breaks = c("Deletion","Duplication","Insertion sequence","None"),name = "Mutation type",na.value = 'white',na.translate = T,guide = guide_legend(order=4,override.aes = list(colour = "black" ))) 

# Evolution experiment themes
parent_colors <- scale_color_manual(values = c(hues::iwanthue(3,hmin=5,hmax=25,plot=F),hues::iwanthue(3,hmin=230,hmax=240,plot=F)),breaks = c("Parent - Plate 1","Parent - Plate 2","Parent - Plate 3","Resistant - Plate 1","Resistant - Plate 2","Resistant - Plate 3"),name = c("Isolate - Replicate"),guide = guide_legend(ncol=2,order=2))

parent_color_group <- scale_fill_manual(values = c("red","purple"),labels = c("Parent","Resistant"),breaks = c("parent","resistant"),name = "Isolate type",guide=guide_legend(nrow=1)) 

status_fill <-  scale_fill_manual(values = c("#FC4E07","pink","blue"),guide = guide_legend(ncol=1,order=1),name = "Status of the OmpK36 porin") 

Limit_scale <- scale_color_manual(breaks = c(T,F),values = c("black","#A9A9A9"),labels = c("Over","Under"),name = "Limit of detection",guide=guide_legend(ncol=1))

# Plasmid categories
categories <- c("#7684C0","maroon","lightgray",'black')
plasmid_type <- c("AA018 with blaKPC", "AA018 with blaKPC + AA552 without blaKPC", "AA552 without blaKPC", "Other blaKPC plasmid")
plasmid_fill <- scale_color_manual(values = categories, breaks = plasmid_type,
                                   labels = gsub(" \\+ ", " +<br>",
                                                 gsub("blaKPC", "<i>bla</i><sub>KPC</sub>", plasmid_type)),
                                   name = "Plasmid configuration",
                                   guide = guide_legend(ncol = 1))

plasmid_type_other <- c("AA018", "AA552", "Other plasmid")
plasmid_fill_other <- scale_fill_manual(values =  c("#7684C0","brown","black"),breaks = plasmid_type_other,name = "Plasmid cluster",guide=guide_legend(ncol=1),na.translate=TRUE,na.value='white') 
plasmid_fill_other_key_border <- scale_fill_manual(values =  c("#7684C0","brown","black"),breaks = plasmid_type_other,name = "Plasmid cluster",guide=guide_legend(ncol=1,override.aes = list(colour = "black" )),na.translate=TRUE,na.value='white') 

# Log 10 frequency scale
log_10_experimental_evolution_scale <-  scale_x_log10(
  limits = c(5e-13, 3e-6),
  breaks = c(1e-12, 1e-10, 1e-8,1e-6),
  labels = c("10<sup>\u221212</sup>", "10<sup>\u221210</sup>", "10<sup>\u22128</sup>", "10<sup>\u22126</sup>")
)

# KPC experimental evolution copy number scale
KPC_copy_number_change_scale <- scale_fill_gradient2(low = "blue", mid = "white" ,high = "red",
                                                     name = expression(log[2]~"fold change"),
                                                     midpoint=0,
                                                     na.value = 'white',
                                                     guide = guide_colorbar(order=3),
                                                     limits=c(-1,3))

# Mobility scale
mobility_scale <- scale_fill_manual(
  values = c("brown", "blue", "#C0C0C0", "transparent"),
  breaks = c("conjugative", "mobilizable", "non-mobilizable", "dummy"),
  labels = c("Conjugative", "Mobilizable", "Non-mobilizable", ""),
  name = "Predicted mobility",
  drop = FALSE,
  guide = guide_legend(nrow = 2, ncol = 3, byrow = TRUE, title.position = "top", label.position = "right",
    override.aes = list(
      fill = c("brown", "blue", "#C0C0C0", "transparent"),
      color = c("transparent", "transparent", "transparent", "transparent")
    )
  )
)
 
# Host range scale
color_range <- hues::iwanthue(n = 4,hmin = 0,hmax = 75,plot = F)
host_range_scale <- scale_fill_manual(values = c("darkgray", color_range ),
                                      breaks = c("Actinomycetota,Bacillota,Pseudomonadota","Gammaproteobacteria","Enterobacterales","Enterobacteriaceae","Klebsiella"),
                                      labels =  c("Actinomycetota, Bacillota, Pseudomonadota","Gammaproteobacteria","Enterobacterales","Enterobacteriaceae","Klebsiella"),
                                      name="Observed host range",
                                      guide = guide_legend(nrow=2, title.position = "top", label.position = "right"))

# Log2 continuous scale for supplemental phylogeny
log2_cont_scale <-  scale_fill_gradient2(
  low = "blue",
  mid = "white",
  high = "red",
  midpoint = 0,
  na.value = "lightgrey",
  name = "log\u2082 copy number",
  guide = guide_colorbar(direction = 'vertical',order=5,
    barwidth = unit(0.5, "cm"),
    barheight = unit(3, "cm"))) 

log2_cont_scale_s6 <- scale_fill_gradient2(
  low = "blue", mid = "white", high = "red", midpoint = 0,
  na.value = "lightgrey",name = "log<sub>2</sub> copy number",
  guide = guide_colorbar(direction = "vertical", order = 5, title.position = "top",
                         barwidth = unit(0.25, "cm"), barheight = unit(1.8, "cm"),
                         theme = theme(legend.title = element_markdown(size = 7, hjust = 0.5))))

# Epidemiology of AA552 plasmid scale
AA552_scale <- scale_fill_manual(values = c("blue",'lightblue','lightgray'), breaks = c("AA552 <i>bla</i><sub>KPC</sub> plasmid","AA552 Non-<i>bla</i><sub>KPC</sub> plasmid","Other"),name="AA552 Status")
AA552_scale_color <- scale_color_manual(values = c("blue",'lightblue','lightgray'), breaks = c("AA552 <i>bla</i><sub>KPC</sub> plasmid","AA552 Non-<i>bla</i><sub>KPC</sub> plasmid","Other"),name="AA552 Status")

# Favorite kable
favorite_kable <- function (x){
  x %>% kable(., format = "html", table.attr = "style='width:100%;'",
              row.names = F) %>% kable_styling(bootstrap_options = c("striped",
                                                                     "hover", "condensed", "responsive"))
}

# Table themes
## Main text
### Figure 3 Table. 
mytheme <- ttheme_minimal(base_family = "Arial",
                          core = list(fg_params = list(hjust = 0, x = 0.01,
                                                       fontsize = 18),
                                      padding = unit(c(5, 2.5), "mm")),
                          colhead = list(fg_params = list(hjust = 0, x = 0.01, fontsize = 18,
                                                          fontface = "bold")),
                          rowhead = list(fg_params = list(hjust = 0, x = 0)))


# Figure themes
## Main text figur
### 1. Figure 1 (Histogram)
figure_1_format <- theme(legend.position = "bottom",
                         axis.line         = element_line(linewidth = 0.3, colour = "black"),
                         axis.ticks        = element_line(linewidth = 0.3, colour = "black"),
                         axis.ticks.length = unit(1, "mm"),
                         axis.text.x  = element_text(size = 8, color = "black", hjust = 0.5),
                         axis.text.y  = element_text(size = 8, color = "black", vjust = 0.5),
                         axis.title   = element_markdown(size = 9.5, color = "black"),
                         legend.text  = element_text(size = 8, color = "black"),
                         legend.title = element_text(size = 9.4, color = "black"),
                         legend.key.size = unit(5, "mm"))

### 2. Figure 2 (Histogram)
figure_2_format <- theme(legend.position = "bottom",
                         axis.text.x  = element_text(size = 8.9, color = "black", hjust = 0.5),
                         axis.text.y  = element_text(size = 8.9, color = "black", vjust = 0.5),
                         axis.title   = element_text(size = 10.4, color = "black"),
                         legend.text  = element_text(size = 8.9, color = "black"),
                         legend.title = element_text(size = 10.4, color = "black"),
                         legend.key.size = unit(5, "mm"),
                         legend.box.spacing   = unit(1, "mm"),
                         legend.margin        = margin(0, 0, 0, 0),
                         legend.key.spacing.y = unit(0.5, "mm"))

### 3. Figure 3 GWAS
#### Table
mytheme_GWAS <- ttheme_minimal(base_family = "Arial",
                               core = list(fg_params = list(hjust = 0, x = 0.05,
                                                            fontsize = 18)),
                               colhead = list(fg_params = list(hjust = 0, x = 0.05, fontsize = 20,
                                                               fontface = "bold")),
                               padding = unit(c(5, 5), "mm"))

#### Descriptive plot theme
figure_3_descriptive_plot_theme <- theme(legend.position = "bottom",
                                         axis.text.x  = element_markdown(size = 5.8, color = "black", hjust = 0.5),
                                         axis.text.y  = element_markdown(size = 5.8, color = "black", vjust = 0.5),
                                         axis.title.y = element_markdown(size = 6.4, colour = "black"),
                                         axis.title.x = element_markdown(size = 6.4, colour = "black"),
                                         legend.title = element_text(size = 7, colour = "black", hjust = 0.5),
                                         legend.text  = element_text(size = 6.4, colour = "black"),
                                         legend.key.size  = unit(0.24, "cm"),
                                         legend.key.width = unit(0.24, "cm"))

#### Tree themes# Overall theme
# Overall theme
consistent_theme_GWAS <- theme(legend.position = 'bottom',
                               legend.direction = "horizontal",
                               legend.title.position = 'top',
                               legend.justification = "center",
                               legend.key = element_rect(colour = 'black', linewidth = 0.2),
                               legend.box.spacing = unit(0.00001, "cm"),
                               legend.key.size  = unit(0.22, "cm"),
                               legend.key.width = unit(0.22, "cm"),
                               legend.key.spacing.x = unit(1, "mm"),
                               legend.key.spacing.y = unit(0.5, "mm"),
                               legend.title = element_text(size = 7, hjust = 0.5),
                               legend.text  = element_text(size = 6, hjust = 0,
                                                           margin = margin(l = 3, unit = "pt")),
                               legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
                               legend.spacing.x = unit(0.05, "cm"))

### 4. Figure 4 Format
#### Overall figure theme
figure_4_format <- theme(legend.position = "bottom",
                         axis.text.x  = element_text(size = 8, color = "black", hjust = 0.5),
                         axis.text.y  = element_text(size = 8, color = "black", vjust = 0.5),
                         axis.title.x = element_markdown(size = 9.5, color = "black"),
                         axis.title.y = element_text(size = 9.5, color = "black"),
                         legend.text  = element_text(size = 8, color = "black"),
                         legend.title = element_text(size = 9.5, color = "black"),
                         plot.title   = element_text(size = 9.5, color = "black"),
                         strip.text   = element_blank(),
                         legend.key.size = unit(5, "mm"),
                         legend.box.spacing   = unit(1, "mm"),
                         legend.margin        = margin(0, 0, 0, 0),
                         legend.key.spacing.y = unit(0.5, "mm"))

#### Theme without markdown
figure_4_format_wo_element_markdown <- theme(legend.position = "bottom",
                                             axis.text.x  = element_text(size = 8, color = "black", hjust = 0.5),
                                             axis.text.y  = element_text(size = 8, color = "black", vjust = 0.5),
                                             axis.title   = element_text(size = 9.5, color = "black"),
                                             legend.text  = element_text(size = 8, color = "black"),
                                             legend.title = element_text(size = 9.5, color = "black"),
                                             plot.title   = element_text(size = 9.5, color = "black"),
                                             legend.key.size = unit(5, "mm"),
                                             legend.box.spacing   = unit(1, "mm"),
                                             legend.margin        = margin(0, 0, 0, 0),
                                             legend.key.spacing.y = unit(0.5, "mm"))
### 5. Figure 5 Format 
## Figure A
figure_5_frequency_format <- theme(legend.position = 'bottom',
                                   legend.title.position = 'top',
                                   axis.text.x  = element_markdown(size = 7.5),
                                   axis.text.y  = element_text(size = 7.5),
                                   axis.title.x = element_text(size = 8.5),
                                   legend.title = element_text(size = 8.5),
                                   legend.text  = element_markdown(size = 7.5),
                                   legend.key.size = unit(3.25, "mm"),
                                   legend.key.justification = "center",
                                   axis.line.y  = element_blank(),
                                   legend.spacing.x = unit(1, "mm"),
                                   legend.key.spacing.x = unit(1, "mm"),
                                   legend.margin = margin(0, 0, 0, 0))

## Panel B 
figure_5_grid <- theme(axis.text    = element_markdown(size = 6.8),
                       legend.title = element_text(size = 7.5),
                       legend.text  = element_text(size = 6.125),
                       legend.key.size = unit(3.25, "mm"),
                       plot.background = element_rect(color = 'white'),
                       legend.key.spacing.y = unit(0.03, 'cm'),
                       axis.text.x = element_markdown(angle = 270, hjust = 0, vjust = 0.5),
                       panel.border = element_blank(),
                       panel.grid = element_blank())

## Supplemental Figures

### 2. Supplemental Figure 2: 
format <-   theme(legend.position = "bottom",
                  axis.text =   element_text(size=18,color="black"),
                  axis.title = element_text(size = 22,color="black"),
                  legend.text =   element_text(size=20,color="black"),
                  legend.title = element_text(size = 22,color="black"),
                  plot.title = element_text(size = 24,color="black")
)



### 3. Supplemental Figure 3 - Phylogeny 
#### Phylogeny (legend styling; matches consistent_theme_GWAS)
consistent_theme_sfigure_3 <- theme(legend.position = 'bottom',
                                    legend.direction = "horizontal",
                                    legend.title.position = 'top',
                                    legend.justification = "center",
                                    legend.key = element_rect(colour = 'black', linewidth = 0.2),
                                    legend.box.spacing = unit(0, "cm"),
                                    legend.key.size  = unit(0.22, "cm"),
                                    legend.key.width = unit(0.22, "cm"),
                                    legend.key.spacing.x = unit(1, "mm"),
                                    legend.key.spacing.y = unit(0.5, "mm"),
                                    legend.title = element_markdown(size = 6.5, hjust = 0.5),
                                    legend.text  = element_markdown(size = 5.5, hjust = 0,
                                                                    margin = margin(l = 3, unit = "pt")),
                                    legend.margin = margin(0, 0, 0, 0),
                                    legend.spacing.x = unit(0.05, "cm"))

#### Descriptive plots (matches figure_3_descriptive_plot_theme)
s_figure_3_descriptive_plot_theme <- theme(legend.position = "bottom",
                                           axis.text.x  = element_markdown(size = 5.8, color = "black", hjust = 0.5),
                                           axis.text.y  = element_markdown(size = 5.8, color = "black", vjust = 0.5),
                                           axis.title.y = element_markdown(size = 6.4, colour = "black"),
                                           axis.title.x = element_markdown(size = 6.4, colour = "black"),
                                           legend.title = element_text(size = 6.5, colour = "black", hjust = 0.5),
                                           legend.text  = element_text(size = 5.5, colour = "black"),
                                           legend.key.size  = unit(0.24, "cm"),
                                           legend.key.width = unit(0.24, "cm"),
                                           panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
                                           axis.ticks   = element_line(colour = "black", linewidth = 0.3)) 

### 4. Supplemental Figure - GWAS hits
#### Descriptive plots
s_figure_4_descriptive_plot_theme <- theme(legend.position = "bottom",
                                           axis.text.x  = element_markdown(size = 5.8, colour = "black", hjust = 0.5),
                                           axis.text.y  = element_markdown(size = 5.8, colour = "black", vjust = 0.5),
                                           axis.title.x = element_markdown(size = 6.4, colour = "black"),
                                           axis.title.y = element_markdown(size = 6.4, colour = "black"),
                                           legend.title = element_text(size = 7, colour = "black", hjust = 0.5),
                                           legend.text  = element_text(size = 6, colour = "black"),
                                           legend.key.size  = unit(0.24, "cm"),
                                           legend.key.width = unit(0.24, "cm"),
                                           panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
                                           axis.ticks   = element_line(colour = "black", linewidth = 0.3))

### Most descriptive plots
s_figure_5_format <- theme(legend.text  = element_markdown(size = 6), 
                           legend.key.size = unit(3, "mm"),
                           axis.title   = element_text(size = 7, colour = "black"),
                           axis.text.x  = element_text(size = 6, colour = "black", angle = 45, hjust = 1),
                           axis.text.y  = element_text(size = 6, colour = "black", vjust = 0.5),
                           strip.text   = element_text(size = 6.5),
                           axis.ticks   = element_line(linewidth = 0.3),
                           legend.box.spacing = unit(0.75, "mm"),
                           legend.margin      = margin(0, 0, 0, 0),
                           legend.title = element_markdown(size = 7, margin = margin(b = 2, unit = "pt")),   # title to keys
                           legend.key.spacing.y = unit(1, "mm")                                      # between key rows
)

#### Panel F and G figure
s_figure_5_FG_format <- theme(legend.text  = element_markdown(size = 6), 
                              legend.key.size = unit(3, "mm"),
                              axis.title   = element_text(size = 7, colour = "black"),
                              axis.text.x  = element_text(size = 6, colour = "black", angle = 45, hjust = 1),
                              axis.text.y  = element_text(size = 6, colour = "black", vjust = 0.5),
                              strip.background = element_rect(color = 'black', linewidth = 0.5),
                              strip.text   = element_markdown(size = 6.5),
                              axis.ticks   = element_line(linewidth = 0.3),
                              legend.box.spacing = unit(0.75, "mm"),
                              legend.margin      = margin(0, 0, 0, 0),
                              legend.title = element_markdown(size = 7, margin = margin(b = 2, unit = "pt")),   # title to keys
                              legend.key.spacing.y = unit(1, "mm")                                            # between key rows
)


### 6. Supplemental Figure 6 
consistent_theme_s_figure_6_wo_md <- theme(legend.position = 'bottom',
                                           legend.title.position = 'top',
                                           legend.key = element_rect(colour = 'black', linewidth = 0.2),
                                           legend.box.spacing = unit(0, "cm"),
                                           legend.key.size  = unit(0.25, "cm"),
                                           legend.key.width = unit(0.25, "cm"),
                                           legend.key.spacing.y = unit(0.5, "mm"),
                                           legend.spacing.x = unit(0.1, "cm"),
                                           legend.title = element_text(size = 7, hjust = 0.5,
                                                                       margin = margin(b = 2, unit = "pt")),
                                           legend.text  = element_text(size = 6, hjust = 0))


### 7. Supplemental Figure 7  
#### Descriptive plot format
s_figure_7_descriptive_plot_format <- theme(legend.position = "bottom",
                                            axis.text.x  = element_markdown(size = 6.5, color = "black"),
                                            axis.text.y  = element_markdown(size = 6.5, color = "black"),
                                            axis.title.y = element_markdown(size = 7.5, color = "black"),
                                            axis.title.x = element_markdown(size = 7.5, color = "black"),
                                            legend.text  = element_text(size = 6, color = "black"),
                                            legend.title = element_text(size = 7, color = "black"),
                                            legend.key.size  = unit(0.25, "cm"),          # CHANGED: 3 mm -> 0.25 cm (match)
                                            legend.key.width = unit(0.25, "cm"),          # CHANGED
                                            legend.key.spacing.y = unit(0.5, "mm"),       # CHANGED: match
                                            plot.title   = element_text(size = 7, color = "black"))

#### GWAS format
consistent_theme_sfigure_7 <- theme(legend.position = 'bottom', legend.direction = "horizontal",
                                    legend.justification = "center", legend.key = element_rect(colour = 'white'),
                                    legend.box.spacing = unit(0, "cm"),
                                    legend.key.size = unit(0.25, "cm"), legend.key.width = unit(0.25, "cm"),
                                    legend.title = element_markdown(size = 7, hjust = 0.5),
                                    legend.text  = element_markdown(size = 6, hjust = 0),
                                    legend.margin = margin(t = -0.19, r = 0, b = 0, l = 0, unit = "cm"),
                                    legend.spacing.x = unit(0.05, "cm"), legend.key.spacing.y = unit(0.5, "mm"))

### Supplemental Figure 8 
s_figure_8_format <- theme(legend.position = "bottom",
                           axis.text.x  = element_markdown(size = 8, color = "black", hjust = 0.5),   # 6 -> 7
                           axis.text.y  = element_markdown(size = 8, color = "black", vjust = 0.5),   # 6 -> 7
                           axis.title   = element_markdown(size = 9, color = "black"),                # 7 -> 8
                           legend.text  = element_markdown(size = 8, color = "black"),                # 6 -> 7
                           legend.title = element_blank(),
                           legend.key.size  = unit(0.5, "cm"),
                           legend.key.width = unit(1, "cm"),                                        # longer line in key
                           legend.margin = margin(t = -3, r = 0, b = 0, l = 0),
                           strip.background = element_rect(fill = "lightgray", color = "black", linewidth = 0.3),
                           strip.text.x = element_markdown(size = 5.25, color = 'black'),                # 5.5 -> 5
                           strip.text.y = element_markdown(size = 5.25, color = 'black'),                # 5.5 -> 5
                           panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
                           axis.ticks   = element_line(linewidth = 0.5),
                           panel.spacing.x = unit(0.2, "cm"),
                           panel.spacing.y = unit(0.2, "cm"),
                           strip.clip = "off")

# Suppemental Figure 9. 
s_figure_9_format <- theme(legend.title = element_markdown(size = 7.5),
                           legend.text  = element_markdown(size = 6.5),
                           legend.key.size = unit(0.25, "cm"),
                           legend.key.spacing.y = unit(0.5, "mm"),
                           axis.text.x  = element_markdown(size = 6.5, colour = "black", hjust = 0.5),
                           axis.text.y  = element_markdown(size = 6.5, colour = "black", vjust = 0.5),
                           axis.title.x = element_markdown(size = 7.5, colour = "black"),
                           axis.title.y = element_markdown(size = 7.5, colour = "black"),
                           axis.ticks   = element_line(linewidth = 0.4),
                           plot.margin  = margin(t = 10, r = 2, b = 1, l = 2, unit = "pt"))   # CHANGED: room above each panel for its letter

s_figure_9_format_plasmid_x <- s_figure_9_format +
  theme(axis.text.x = element_markdown(size = 6, colour = "black", angle = 45, hjust = 1, vjust = 1))


### 10. Supplemental Figure 10 
s_figure_10_format <- theme(legend.position = "bottom",
                            axis.text  = element_text(size = 7, color = "black"),
                            axis.title = element_text(size = 8, color = "black"),
                            plot.title = element_text(size = 8, color = "black"),
                            legend.title = element_text(size = 10, color = "black"),
                            legend.text  = element_text(size = 9, color = "black"),
                            legend.key.size = unit(0.45, "cm"),
                            legend.key.spacing.y = unit(1, "mm"),
                            legend.margin = margin(t = 0, b = 0, unit = "cm"),
                            plot.margin  = margin(t = 10, r = 2, b = 1, l = 2, unit = "pt"))   # CHANGED: room above each panel for its letter

### 11. Supplemental Figure 11 
s_figure_11_format_A <- theme(legend.title.position = 'top',
                              legend.title = element_text(size = 9, colour = "black"),
                              legend.text  = element_text(size = 8, colour = "black"),
                              legend.key.size = unit(0.4, "cm"),
                              legend.key.spacing.y = unit(0.5, "mm"),
                              axis.text.x  = element_text(size = 8, colour = "black"),
                              axis.text.y  = element_text(size = 8, colour = "black"),
                              axis.title   = element_text(size = 9, colour = "black"),
                              strip.text   = element_text(size = 8),
                              axis.ticks   = element_line(linewidth = 0.3),
                              plot.margin  = margin(t = 12, r = 3, b = 2, l = 3, unit = "pt"))

s_figure_11_format_BCDE <- theme(axis.text.x  = element_markdown(angle = 45, hjust = 1, size = 8, colour = "black"),
                                 axis.text.y  = element_markdown(size = 8, colour = "black"),
                                 axis.title.y = element_markdown(size = 9, colour = "black"),
                                 strip.text   = element_text(size = 7),
                                 legend.title.position = 'top',
                                 legend.title = element_text(size = 9, colour = "black"),
                                 legend.text  = element_text(size = 8, colour = "black"),
                                 legend.key.size = unit(0.4, "cm"),
                                 axis.ticks   = element_line(linewidth = 0.3),
                                 plot.margin  = margin(t = 12, r = 3, b = 2, l = 3, unit = "pt"))

### 12. Supplementary Figure 12
s_figure_12_format_B <- theme(legend.position = "bottom", legend.box = "horizontal",
                              legend.title = element_text(size = 7.5, colour = "black"),
                              legend.text = element_markdown(size = 6.5, colour = "black"),
                              legend.key.size = unit(0.35, "cm"), legend.key.spacing.y = unit(0.4, "mm"),
                              plot.margin = margin(t = 12, r = 3, b = 2, l = 3, unit = "pt"))

### 13. Supplementary Figure 13 
s_figure_13_format <- theme(axis.text.x  = element_text(size = 6.5, colour = "black"),
                            axis.text.y  = element_markdown(size = 6.5, colour = "black"),
                            axis.title   = element_text(size = 7.5, colour = "black"),
                            strip.text   = element_text(size = 6.5),
                            axis.ticks   = element_line(linewidth = 0.3),
                            legend.position = "none",
                            plot.margin  = margin(t = 12, r = 3, b = 2, l = 3, unit = "pt"))

### 14. Supplemental Figure 14 - Genotype count 
s_figure_14_format <-  theme(legend.position = "bottom",
                             axis.text =   element_text(size=10,color="black"),
                             axis.title = element_text(size = 12,color="black"),
                             legend.text =   element_text(size=10,color="black"),
                             legend.title = element_text(size = 12,color="black"),
                             plot.title = element_text(size = 12,color="black"),legend.key.size = unit(0.4, "cm")
                             )
