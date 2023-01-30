require(biomaRt)
require(SingleCellExperiment)
require(liana)
require(tibble)
require(purrr)
require(Seurat)
require(reticulate)
require(tidyverse)
require(CrossTalkeR)

library(tidyverse)
library(magrittr)
library(liana)
library(CrossTalkeR)
library(dplyr)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(zellkonverter)
library(nichenetr)
library(tidyverse)
library(DelayedArray)


###EXTRACTING DATA
require(tibble)
require(biomaRt)
require(tidyr)
require(dplyr)
require(SingleCellExperiment)

##############################################################  MAPPING RAT TO HUMAN 
library(tidyverse)
library(OmnipathR)
library(liana)
library(magrittr)

show_homologene()

data1 <- readRDS("scvi_output_clustered_subs.rds")

data = as.SingleCellExperiment(data1)


case <-assays(data)[["logcounts"]][,colData(data)$cc=='case']

control <-assays(data)[["logcounts"]][,colData(data)$cc!='case']

metadatacase <-tibble::tibble(rownames(colData(data))[colData(data)$cc=='case'],
                              as.character(colData(data)$leiden0.3[colData(data)$cc=='case']))
metadatacontrol<- tibble::tibble(rownames(colData(data))[colData(data)$cc!='case'],
                                 as.character(colData(data)$leiden0.3[colData(data)$cc!='case']))


human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
rat <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
#listAttributes(rat)

#exclude 
#genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = rownames(case) , mart = rat, attributesL = c("hgnc_symbol","hgnc_id",'ensembl_gene_id'), martL = human, uniqueRows=T)

genesV2 = getLDS(attributes = c("rgd_symbol"), 
                 filters = "rgd_symbol", 
                 values = rownames(case) , 
                 mart = rat, 
                 attributesL = c("hgnc_symbol"), 
                 martL = human, 
                 uniqueRows=T)


mapping <- genesV2[genesV2$HGNC.symbol!='',]

mappinglist <- mapping$HGNC.symbol
names(mappinglist) <- mapping$RGD.symbol

rownames(case) <- mappinglist[rownames(case)]
case <- case[!is.na(rownames(case)),]
rownames(control) <- mappinglist[rownames(control)]
control <- control[!is.na(rownames(control)),]
dim(case)
dim(control)

##CASE SEURAT OBJECT FOR LIANA

lia_case = CreateSeuratObject(counts = case)
Idents(lia_case) = metadatacase$`as.character(colData(data)$leiden0.3[colData(data)$cc == "case"])`#metadatacase$`as.character(...)`
head(lia_case@meta.data)



#CONTROL SEURAT OBJECT FOR LIANA

lia_cont = CreateSeuratObject(counts = control)
Idents(lia_cont) = metadatacontrol$`as.character(colData(data)$leiden0.3[colData(data)$cc != "case"])`#metadatacontrol$`as.character(...)`
head(lia_cont@meta.data)


#RUNNING LIANA FOR CASE AND CONTROL
#Idents(lia_case) = lia_case@meta.data$idents
lia_case = NormalizeData(lia_case, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test <-  liana_wrap(lia_case)



lia_cont = NormalizeData(lia_cont, normalization.method = "LogNormalize", scale.factor = 10000)
liana_test2 <- liana_wrap(lia_cont)

paths1 <-liana_test$cellphonedb%>%
  filter(pvalue<=0.01) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 

paths2 <-liana_test2$cellphonedb %>% 
  filter(pvalue<=0.01) %>%
  mutate(gene_A=gsub('_','*',ligand.complex)) %>%
  mutate(gene_B=gsub('_','*',receptor.complex)) %>%
  mutate(source=gsub("_","",source)) %>%
  mutate(target=gsub("_","",target)) %>%
  mutate(type_gene_A = rep('Ligand',length(.data$ligand))) %>%
  mutate(type_gene_B = rep('Receptor',length(.data$receptor))) 

teste <- list()
teste[['control']] <- paths2 
teste[['case']] <- paths1

teste
#
library(rmarkdown)
library(tinytex)

data <- generate_report(lrpaths = teste,
                        genes=NULL,
                        out_path='/Users/newuser/sidrah/stent/',
                        threshold=0,
                        out_file = 'report_rat_liana.pdf',
                        output_fmt = "pdf_document",
                        report = T,
                        sel_columns = c('source','target','gene_A','gene_B','type_gene_A','type_gene_B','lr.mean'))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++ SANKEY PLOTS ++++++++++++++++++++++++++++++++

all_data <- readRDS("LR_data_final.Rds")


for(i in 1:length(names(all_data@pca))){
  curr <- names(all_data@pca)[i]
  if(str_detect(curr, '_x_', negate = FALSE) & str_detect(curr, '_ggi', negate = FALSE)){
    rmd_title <- paste0(curr,'_tbl')
    rmd_title1 <- paste0(curr,'_pca')
    x <- max(abs(all_data@pca[[curr]]$x[,1]))
    y <- max(abs(all_data@pca[[curr]]$x[,2]))
    z_x <- all_data@pca[[curr]]$x[,1]
    z_y <- all_data@pca[[curr]]$x[,2]
    ver_zx <- ifelse(abs(z_x)>=(2*all_data@pca[[curr]]$sdev[1]),1,0)
    ver_zy <- ifelse(abs(z_y)>=(2*all_data@pca[[curr]]$sdev[2]),1,0)
    print(fviz_pca_biplot(all_data@pca[[curr]],
                          axes = c(1,2),
                          pointshape = 21, pointsize = 0.5,labelsize = 10,
                          repel = FALSE,max.overlaps=100,label='var')+
            geom_label_repel(aes(label=ifelse((ver_zx | ver_zy),rownames(all_data@pca[[curr]]$x),NA)),size = 5)+
            xlim(-x, x)+
            ylim(-y, y)+
            ggtitle(curr)+
            theme(text = element_text(size = 7.5),
                  axis.title = element_text(size = 7.5),
                  axis.text = element_text(size = 7.5)))
  }  
}




for(i in 1:length(names(all_data@pca))){
  curr <- names(all_data@pca)[i]
  if(str_detect(curr, '_x_', negate = FALSE) & str_detect(curr, '_ggi', negate = FALSE)){
    rmd_title <- paste0(curr,'_tbl')
    rmd_title1 <- paste0(curr,'_pca')
    x <- max(abs(all_data@pca[[curr]]$x[,1]))
    y <- max(abs(all_data@pca[[curr]]$x[,2]))
    z_x <- all_data@pca[[curr]]$x[,1]
    z_y <- all_data@pca[[curr]]$x[,2]
    ver_zx <- ifelse(abs(z_x)>=(2*all_data@pca[[curr]]$sdev[1]),1,0)
    ver_zy <- ifelse(abs(z_y)>=(2*all_data@pca[[curr]]$sdev[2]),1,0)
    print(fviz_pca_biplot(all_data@pca[[curr]],
                          axes = c(1,2),
                          pointshape = 21, pointsize = 0.5,labelsize = 10,
                          repel = FALSE,max.overlaps=100,label='var')+
            geom_label_repel(aes(label=ifelse((ver_zx | ver_zy),rownames(all_data@pca[[curr]]$x),NA)),size = 5)+
            xlim(-x, x)+
            ylim(-y, y)+
            ggtitle(curr)+
            theme(text = element_text(size = 7.5),
                  axis.title = element_text(size = 7.5),
                  axis.text = element_text(size = 7.5)))
  }  
} 





options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='EPHB4',
            threshold = 20)

ggplot2::ggsave("EPHB4_sankey_targetted_rat.pdf",width = 15) 

#### Plotting all collagen interactions
res2 <- readRDS("LR_data_final.Rds")

l = unique(grep("COL", res2@tables$case_x_control$Ligand, value = TRUE))

#1 COL12A1
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL12A1',
            threshold = 20)

#2
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL14A1',
            threshold = 20)

#3
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL15A1',
            threshold = 20)

#4
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL16A1',
            threshold = 20)

#5
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL1A1',
            threshold = 20)

#6
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL3A1',
            threshold = 20)

#7
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL4A1',
            threshold = 20)

#8
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL4A2',
            threshold = 20)

#9
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL4A3',
            threshold = 20)

#10
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL4A6',
            threshold = 20)


#11
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL5A1',
            threshold = 20)

#12
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL5A2',
            threshold = 20)

#13
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL6A1',
            threshold = 20)

#14
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL6A2',
            threshold = 20)

#15
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL6A5',
            threshold = 20)

#16
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL6A6',
            threshold = 20)

#17
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL8A1',
            threshold = 20)

#18
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL8A2',
            threshold = 20)

#19
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL18A1',
            threshold = 20)

#20
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL5A3',
            threshold = 20)

#21
options(repr.plot.width=15)
plot_sankey(res2@tables$case_x_control,
            target='COL13A1',
            threshold = 20)



##++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




require(CrossTalkeR)
require(stringr)
require(ggplot2)
require(dplyr)
require(colorBlindness)
require(patchwork)
require(igraph)

pdf('CCIs.pdf',width = 24,height=8)
par(mfrow=c(1,3))
plot_cci(data@graphs$control,
         coords = data@coords[V(data@graphs$control)$name,],
         colors = data@colors[V(data@graphs$control)$name],
         plt_name = "control",
         pg=TRUE,
         leg=TRUE,
         emax=max(E(data@graphs$control)$LRScore))
plot_cci(data@graphs$case,
         coords = data@coords[V(data@graphs$case)$name,],
         colors = data@colors[V(data@graphs$case)$name],
         plt_name = "case",
         pg=TRUE,
         leg=TRUE,
         emax=max(E(data@graphs$control)$LRScore))
plot_cci(data@graphs$case_x_control,
         coords = data@coords[V(data@graphs$case_x_control)$name,],
         colors = data@colors[V(data@graphs$case_x_control)$name],
         plt_name = "case_x_control",
         pg=TRUE,
         leg=TRUE,emax=max(E(data@graphs$control)$LRScore))
dev.off()

####++++++++++++++++++++++++++++++++++++++PAGERANK CONTROL AND CASE
####+
####+
#+++++++CONTROL
library(forcats)
library(dplyr)

cont = data@rankings$control

cont %>%
  mutate(name = reorder(nodes, Pagerank)) %>%
  ggplot( aes(x=Pagerank, y=fct_reorder(nodes, Pagerank))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  xlab("Pagerank") + ggtitle('Pagerank log ratio for control CCI') + ylab('reorder(nodes, Pagerank)')+
  theme_bw()

#### CASE 

cas_page = data@rankings$case

cas_page %>%
  mutate(name = reorder(nodes, Pagerank)) %>%
  ggplot( aes(x=Pagerank, y=fct_reorder(nodes, Pagerank))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  xlab("Pagerank") + ggtitle('Pagerank log ratio for case CCI') + ylab('reorder(nodes, Pagerank)')+
  theme_bw()



###++++++ CGI

cont_gg = data@rankings$control_ggi

cont_gg %>%
  mutate(name = reorder(nodes, Pagerank)) %>%
  top_n(n = 10, wt = Pagerank) %>%
  ggplot( aes(x=Pagerank, y=fct_reorder(nodes, Pagerank))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  xlab("Pagerank") + ggtitle('Pagerank log ratio for control GGI') + ylab('reorder(nodes, Pagerank)')+
  theme_bw()

#### CASE 

cas_gg = data@rankings$case_ggi

cas_gg %>%
  mutate(name = reorder(nodes, Pagerank)) %>%
  top_n(n = 10, wt = Pagerank) %>%
  ggplot( aes(x=Pagerank, y=fct_reorder(nodes, Pagerank))) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  xlab("Pagerank") + ggtitle('Pagerank log ratio for case GGI') + ylab('reorder(nodes, Pagerank)')+
  theme_bw()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++GSEA 

library(clusterProfiler)
complrall <- res2@tables$case_x_control %>%
  dplyr::group_by(Ligand.Cluster,Receptor.Cluster) %>%
  summarise(gene_A =list(Ligand),gene_B = list(Receptor),LRscore = list(LRScore))

complr <- complrall[grepl('Fib',complrall$Ligand.Cluster)|
                      grepl('Fib',complrall$Receptor.Cluster),]

print(complr, n = 30)

aux <- tibble(Description='',ratio=0,cci='')

for(i in 1:dim(complr)[1]){
  selup <- complr[i,]$LRscore[[1]] < 0
  if(sum(selup) != 0){
    neggenes <- union(complr[i,]$gene_A[[1]][selup],complr[i,]$gene_B[[1]][selup])
    neggenes <- bitr(neggenes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
    resdown <- enrichKEGG(neggenes$ENTREZID)
    dn <- resdown@result %>%
      filter(p.adjust<=0.05)
  }else{
    dn = NULL
  }
  selup <- complr[i,]$LRscore[[1]] > 0
  if(sum(selup) != 0){
    posgenes <- union(complr[i,]$gene_A[[1]][selup],complr[i,]$gene_B[[1]][selup])
    posgenes <- bitr(posgenes,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = "org.Hs.eg.db")
    resup <- enrichKEGG(posgenes$ENTREZID)
    up <- resup@result %>%
      filter(p.adjust<=0.05)
  }else{
    up = NULL
  }
  if(!is.null(dn) & !is.null(up)){
    joined <- merge(x=up,y=dn,by.x="Description",by.y="Description",all=TRUE)
    joined$pvalue.x[is.na(joined$pvalue.x)]=1
    joined$pvalue.y[is.na(joined$pvalue.y)]=1
    joined$ratio <- -log10(joined$p.adjust.x/joined$p.adjust.y)
    joint <- joined[,c("Description","ratio")]
    joint$cci = paste0(complr[i,]$Ligand.Cluster,"|",complr[i,]$Receptor.Cluster)
  }else if(is.null(dn) & !is.null(up)){
    joined <-  up
    joined$ratio <- -log10(joined$p.adjust/1)
    joint <- joined[,c("Description","ratio")]
    joint$cci = paste0(complr[i,]$Ligand.Cluster,"|",complr[i,]$Receptor.Cluster)
  } else if(!is.null(dn) & !is.null(up)){
    joined <-  down
    joined$ratio <- -log10(1/joined$p.adjust)
    joint <- joined[,c("Description","ratio")]
    joint$cci = paste0(complr[i,]$Ligand.Cluster,"|",complr[i,]$Receptor.Cluster)
  }
  aux <- bind_rows(aux,joint)
}


final <- aux[2:dim(aux)[1],] %>%
  reshape2::dcast(Description~cci,value.var = "ratio",fun.aggregate = mean,)
final[is.na(final)] <-0
tmp = final[,2:dim(final)[2]]
rownames(tmp) <- final[,1]
seledg <- sub('_','|',data@stats$case_x_control[data@stats$case_x_control$p<0.05,]$columns_name)
tmp <- tmp[apply(tmp, 1,var) > summary(apply(tmp, 1,var))[[5]]
           ,match(seledg,colnames(tmp),nomatch = F)]

pdf(file = "ec2.pdf",   # The directory you want to save the file in
    width = 25, # The width of the plot in inches
    height = 35)

library(circlize)
f1 = colorRamp2(seq(-max(abs(tmp)), max(abs(tmp)), length = 3), c("blue", "#EEEEEE", "red"))
rownames(tmp) <- str_wrap(rownames(tmp), width = 20)
colnames(tmp) <- str_wrap(colnames(tmp), width = 17)
ComplexHeatmap::Heatmap(tmp,
                        col=f1,
                        name='-log(pvalue_up/pvalue_down)',clustering_method_column='average'
                        ,clustering_method_row='average', column_names_gp = grid::gpar(fontsize = 25),
                        row_names_gp = grid::gpar(fontsize = 20)) 
dev.off()


######## MERGED

tmp_ec2 = tmp
tmp_ec1 = tmp
tmp_fib = tmp
tmp_vsmc = tmp



###VSMC

l1 = c('AGE-RAGE signaling\npathway in diabetic\ncomplications', 'Cholesterol\nmetabolism', 'Notch signaling\npathway', 
       'Fluid shear stress\nand atherosclerosis', 'Cytokine-cytokine\nreceptor interaction',
       'Regulation of actin\ncytoskeleton', 'TGF-beta signaling\npathway', 'MAPK signaling\npathway', 'Cell adhesion\nmolecules',
       'Leukocyte\ntransendothelial\nmigration', 'Rap1 signaling\npathway' , 'PI3K-Akt signaling\npathway', 'Focal adhesion','ECM-receptor\ninteraction')

vsmc = tmp_vsmc[l1,]

vsmc['pathway'] = rownames(vsmc)
vsmc

### FIBROBLAST

rownames(tmp_fib)

l2 = c( 'Cholesterol\nmetabolism', 'Ras signaling\npathway', 'TGF-beta signaling\npathway', 
        'Hippo signaling\npathway', 'MAPK signaling\npathway', 'AGE-RAGE signaling\npathway in diabetic\ncomplications', 'Rap1 signaling\npathway',
        'Fluid shear stress\nand atherosclerosis', 'Cytokine-cytokine\nreceptor interaction', 'Cell adhesion\nmolecules', 'Regulation of actin\ncytoskeleton',
        'Focal adhesion', 'PI3K-Akt signaling\npathway', 'ECM-receptor\ninteraction')


fib = tmp_fib[l2,]
fib['pathway'] = rownames(fib)

fib
###EC1

rownames(tmp_ec1)
l3 = c('Notch signaling\npathway', 'Cholesterol\nmetabolism', 'Ras signaling\npathway', 'AGE-RAGE signaling\npathway in diabetic\ncomplications',
       'Cell adhesion\nmolecules', 'MAPK signaling\npathway', 'Rap1 signaling\npathway', 'Fluid shear stress\nand atherosclerosis', 
       'Focal adhesion', 'PI3K-Akt signaling\npathway', 'ECM-receptor\ninteraction')


ec1 = tmp_ec1[l3,]
ec1['pathway'] = rownames(ec1)
ec1
##EC2

rownames(tmp_ec2)

l4 = c('ECM-receptor\ninteraction', 'AGE-RAGE signaling\npathway in diabetic\ncomplications', 'Leukocyte\ntransendothelial\nmigration', 
       'Cell adhesion\nmolecules', 'Th1 and Th2 cell\ndifferentiation', 'Notch signaling\npathway', 'Fluid shear stress\nand atherosclerosis',
       'Hippo signaling\npathway', 'Ras signaling\npathway', 'Rap1 signaling\npathway', 'PI3K-Akt signaling\npathway', 'Focal adhesion', 
       'Regulation of actin\ncytoskeleton')


ec2 = tmp_ec2[l4,]
ec2['pathway'] = rownames(ec2)
ec2

df1 = merge(x = vsmc, y = fib, by = 'pathway',all = TRUE)
df2 = merge(x = df1, y = ec1, by = 'pathway', all = TRUE)
df_final = merge(x= df2, y = ec2, by = 'pathway', all = TRUE)
df_final

rownames(df_final) = df_final$pathway

df_final2 = subset(df_final, select = -c(pathway))
df_final2[is.na(df_final2)] <- 0
colnames(df_final2)
df_final2


dff_final = df_final2[, c('3 : Fibroblast|\n5 : VSMC.x', '3 : Fibroblast|\n0 : Bcells1',  '3 : Fibroblast|\n6 : Macrophage1',
                          '3 : Fibroblast|\n10 : Macrophage2', '3 : Fibroblast|\n4 : Bcells2', '3 : Fibroblast|\n7 : EC1.x', '3 : Fibroblast|\n2 : NKcell1',
                          '3 : Fibroblast|\n5 : VSMC.y','3 : Fibroblast|\n7 : EC1.y', '7 : EC1|3 :\nFibroblast.x', '9 : EC2|3 :\nFibroblast.x',
                          '9 : EC2|3 :\nFibroblast.y', '4 : Bcells2|3 :\nFibroblast', '7 : EC1|3 :\nFibroblast.y', '3 : Fibroblast|\n3 : Fibroblast',
                          
                          '5 : VSMC|2 :\nNKcell1',
                          '5 : VSMC|6 :\nMacrophage1', '5 : VSMC|7 : EC1.y', '5 : VSMC|4 :\nBcells2', '5 : VSMC|7 : EC1.x',
                          '5 : VSMC|10 :\nMacrophage2', '4 : Bcells2|5 :\nVSMC','9 : EC2|5 : VSMC.x', '5 : VSMC|5 : VSMC',
                          
                          '7 : EC1|10 :\nMacrophage2', '7 : EC1|6 :\nMacrophage1', '7 : EC1|7 : EC1',
                          
                          '9 : EC2|0 :\nBcells1', '9 : EC2|11 :\nMacrophage3', '9 : EC2|1 : Tcell',
                          '9 : EC2|8 :\nNKcell2', '9 : EC2|10 :\nMacrophage2',
                          '9 : EC2|4 :\nBcells2', '9 : EC2|5 : VSMC.y',  '4 : Bcells2|9 :\nEC2','9 : EC2|9 : EC2')]


colnames(dff_final) = c('Fib| VSMC.x', 'Fib| Bcells1',  'Fib| Macro1',
                        'Fib| Macro2', 'Fib| Bcells2', 'Fib| EC1.x', ' Fib| NKcell1',
                        'Fib| VSMC.y','Fib| EC1.y', 'EC1| Fib.x', 'EC2| Fib.x',
                        'EC2| Fib.y', 'Bcells2| Fib', 'EC1| Fibt.y','Fib| Fib',
                        
                        'VSMC| NKcell1',
                        'VSMC| Macro1', 'VSMC| EC1.y', 'VSMC| Bcells2', 'VSMC| EC1.x',
                        'VSMC| Macro2', 'Bcells2| VSMC','EC2| VSMC.x', 'VSMC| VSMC',
                        
                        'EC1| Macro2',  'EC1| Macro1','EC1| EC1',
                        
                        'EC2| Bcells1', 'EC2| Macro3', 'EC2| Tcell',
                        'EC2| NKcell2', 'EC2| Macro2',
                        'EC2| Bcells2', 'EC2| VSMC.y',  'Bcells2| EC2','EC2| EC2')



pdf(file = "GSEA.pdf",   # The directory you want to save the file in
    width = 35, # The width of the plot in inches
    height = 30)

library(circlize)
f1 = colorRamp2(seq(-max(abs(dff_final)), max(abs(dff_final)), length = 3), c("blue", "#EEEEEE", "red"))
rownames(dff_final) <- str_wrap(rownames(dff_final), width = 20)
colnames(dff_final) <- str_wrap(colnames(dff_final), width = 17)
ComplexHeatmap::Heatmap(dff_final,
                        col=f1,cluster_rows = FALSE, cluster_columns = FALSE,
                        name='-log(pvalue_up/pvalue_down)', column_names_gp = grid::gpar(fontsize = 25),
                        row_names_gp = grid::gpar(fontsize = 20))

dev.off()






##### top interactions
####unique interactions


all <- res2@tables$case_x_control %>%
  dplyr::group_by(Ligand.Cluster,Receptor.Cluster)

complr_fib <- all[grepl('Fib',all$Receptor.Cluster),]

complr_fib$Ligand[unique(complr_fib$Ligand)]




complr_ec1 <- complr_fib[grepl('EC1',complr_fib$Ligand.Cluster),]
complr_ec2 <- complr_fib[grepl('EC2',complr_fib$Ligand.Cluster),]
complr_vsmc <- complr_fib[grepl('VSMC',complr_fib$Ligand.Cluster),]
complr_bcell = complr_fib[grepl('Bcells1',complr_fib$Ligand.Cluster),]


complr_tcell <- complr_fib[grepl('Tcell',complr_fib$Ligand.Cluster),]
complr_macro2 <- complr_fib[grepl('Macrophage2',complr_fib$Ligand.Cluster),]
complr_macro3 <- complr_fib[grepl('Macrophage3',complr_fib$Ligand.Cluster),]
complr_nk = complr_fib[grepl('NKcell1',complr_fib$Ligand.Cluster),]


complr_fibro <- complr_fib[grepl('Fibroblast',complr_fib$Ligand.Cluster),]
complr_bcell2 <- complr_fib[grepl('Bcells2',complr_fib$Ligand.Cluster),]
complr_nk2 <- complr_fib[grepl('NKcell2',complr_fib$Ligand.Cluster),]
complr_macr1 = complr_fib[grepl('Macrophage1',complr_fib$Ligand.Cluster),]

# saving interactions

write.csv(complr_ec1, file = 'ec1_intraction.csv')
write.csv(complr_ec2, file = 'ec2_interaction.csv')
write.csv(complr_vsmc, file = 'vsmc_interaction.csv')

write.csv(complr_bcell, file = 'bcell_intraction.csv')
write.csv(complr_bcell2, file = 'bcell2_interaction.csv')
write.csv(complr_fibro, file = 'fibro_interaction.csv')

write.csv(complr_macr1, file = 'macro1_intraction.csv')
write.csv(complr_macro2, file = 'macro2_interaction.csv')
write.csv(complr_macro3, file = 'macro3_interaction.csv')

write.csv(complr_nk, file = 'nk_intraction.csv')
write.csv(complr_nk2, file = 'nk2_interaction.csv')
write.csv(complr_tcell, file = 'tcell_interaction.csv')

write.csv(complr_fib, file = 'fibro_all.csv')


#########+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
