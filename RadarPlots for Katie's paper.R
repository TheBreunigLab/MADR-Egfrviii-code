library(Matrix)
library(Seurat)
library("anndata")
library(harmony)
library(dplyr)
library(ggplot2)
library(ggradar)

setwd('/media/david/4TBNvMe/scRNAseq/Analysis/Adult Glioma/Katie hexagonal plot/Human')
data <- read_h5ad("/media/david/4TBNvMe/scRNAseq/Analyzed datasets/Human EgfrvIII filtered subtypes with clusters named 111925.h5ad")
Adult <- CreateSeuratObject(counts = t(data$layers[["counts"]]), meta.data = data$obs)
Adult <- NormalizeData(Adult)
Adult <- FindVariableFeatures(Adult)
Adult <- ScaleData(Adult,  features=rownames(Adult))
Adult <- RunPCA(Adult, npcs=100)
Adult <- RunHarmony(Adult, group.by.vars = 'sample_label', max.iter.harmony=50)
Adult <- FindNeighbors(Adult,dims = 1:45)
Adult <- RunUMAP(Adult, reduction = "harmony", dims = 1:55)
Adult <- FindClusters(Adult, resolution=0.3, algorithm = 4,method = "igraph")
p1<- DimPlot(Adult, raster=FALSE,group.by='Final_CellType', label=T)+ggtitle('Leiden')
p2<- DimPlot(Adult, group.by = 'subtype', raster=FALSE)
p3<- FeaturePlot(Adult, features = c('EGFR','TOP2A', 'OLIG2','CSF1R', 'PECAM1', 'ALDH1L1', 'MOG', 'GFAP', 'DCX', 'CD8A', 'CD4'),
                 ncol = 4, order = T, cols = c('#fbf4ec', 'red', '#b10000'),raster=FALSE)
(p2|p1)/p3

Original <- Adult
p2<- DimPlot(Adult, group.by = 'subtype', raster=FALSE, cols = c('#1fe98b','#5ba1d3','#fd9f4b'))
(p2|p1)/p3


Adult <- subset(Adult, subset = Final_CellType==c("Tumor Astro-like","Tumor Cycling", "Tumor Mesenchymal-like",
                                                  "Tumor Neuron-like", "Tumor OPC-like 1", "Tumor OPC-like 2", 
                                                  "Tumor OPC-like 3"))
p4<- DimPlot(Adult, raster=FALSE,group.by='Final_CellType')+ggtitle('Leiden')

#################### Mouse
data <- read_h5ad('/media/david/4TBNvMe/scRNAseq/Analyzed datasets/All Katie samples together clean non-processed more samples.h5ad')

MADR <- CreateSeuratObject(counts = t(data$X), meta.data = data$obs)
MADR <- NormalizeData(MADR)
MADR <- FindVariableFeatures(MADR)
MADR <- ScaleData(MADR)
MADR <- RunPCA(MADR, npcs=100)
MADR <- RunHarmony(MADR, group.by.vars = 'sample_label', max.iter.harmony=50)
MADR <- FindNeighbors(MADR)
MADR <- RunUMAP(MADR, reduction = "harmony", dims = 1:75)
MADR <- FindClusters(MADR, resolution=2)
MADR <- subset(MADR, subset = sample_label==c('Egfr glioma 1 hit','Egfr glioma 3 hit female 1', 'Egfr glioma 3 hit male', 'Egfr glioma 4 hit Nf1 tissue'))
p1<- DimPlot(MADR, label = T)+NoLegend()+ggtitle('Leiden')
p2<- DimPlot(MADR, group.by = 'sample_label')
p3<- FeaturePlot(MADR, features = c('Top2a', 'Mog','Gfap','Pdgfra','Olig2','Cd3g','Csf1r', 'Cldn5', 'Pdgfrb'),
                 ncol = 5, order = T, cols = c('lightblue', 'lightgoldenrod1', 'coral1'))
(p2|p1)/p3

Original_m <- MADR
MADR <- subset(MADR, idents = c('25','26','27','21','7', '28'), invert=T)

p4<- DimPlot(MADR, label = T)+NoLegend()+ggtitle('Filtered')
p5<- DimPlot(Original_m, label = T)+NoLegend()+ggtitle('Original') 
p5|p4
rm(p1,p2,p3,p4,p5)
##-#####################-##
#### Add Module Scores ####
##-#####################-##
# Human
 # NPC1 Module Score ####
cd_features <- list(c("DLL3",'DLL1','SOX4','TUBB3','HES6','TAGLN3','NEU4','MARCKSL1','CD24','STMN1','TCF12','BEX1','OLIG1','MAP2','FXYD6','PTPRS','MLLT11','NPPA','BCAN','MEST','ASCL1','BTG2','DCX','NXPH1','HN1','PFN2','SCG3','MYT1','CHD7','ADGRG1','TUBA1A','PCBP4','ETV1','SHD','TNR','AMOTL2','DBN1','HIP1','ABAT','ELAVL4','LMF1','GRIK2','SERINC5','TSPAN13','ELMO1','GLCCI1','SEZ6L','LRRN1','SEZ6','SOX11'))
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'NPC', search = F)
 # MES1 module score ####
cd_features <- list(c('CHI3L1', 'ANXA2', 'ANXA1', 'CD44', 'VIM', 'MT2A', 'C1S', 'NAMPT', 'EFEMP1', 'C1R', 'SOD2', 'IFITM3', 'TIMP1', 'SPP1', 'A2M', 'S100A11','MT1A', 'S100A10', 'FN1', 'LGALS1', 'S100A16', 'CLIC1', 'MGST1', 'RCAN1', 'TAGLN2', 'NPC2', 'SERPING1', 'C8orf4', 'EMP1', 'APOE', 'CTSB', 'C3', 'LGALS3','EMP3', 'SERPINA3', 'ACTN1', 'PRDX6', 'IGFBP7', 'SERPINE1', 'PLP2', 'MGP', 'CLIC4', 'GFPT2', 'GSN', 'NNMT', 'TUBA1C', 'GJA1', 'TNFRSF1A', 'WWTR1'))
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'MES1', search = F)
 # MES2 module score ####
cd_features <- list(c('HILPDA', 'ADM', 'DDIT3', 'NDRG1', 'HERPUD1', 'DNAJB9', 'TRIB3', 'ENO2', 'AKAP12', 'SQSTM1', 'MT1A', 'ATF3', 'NAMPT','NRN1', 'SLC2A1','BNIP3', 'LGALS3', 'INSIG2', 'IGFBP3', 'PPP1R15A', 'VIM', 'PLOD2', 'GBE1', 'SLC2A3', 'FTL', 'WARS','ERO1A', 'XPOT', 'HSPA5', 'GDF15', 'ANXA2', 'EPAS1', 'LDHA', 'P4HA1', 'SERTAD1', 'PFKP', 'PGK1', 'EGLN3', 'SLC6A6', 'CA9', 'BNIP3L', 'RPL21', 'TRAM1', 'UFM1', 'ASNS', 'GOLT1B', 'ANGPTL4', 'SLC39A14', 'CDKN1A', 'HSPA9'))
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'MES2', search = F)
 # MES1+2 module score ####
cd_features <- list(c('CHI3L1', 'ANXA2', 'ANXA1', 'CD44', 'VIM', 'MT2A', 'C1S', 'NAMPT', 'EFEMP1', 'C1R', 'SOD2', 'IFITM3', 'TIMP1', 'SPP1','A2M', 'S100A11','MT1A', 'S100A10', 'FN1', 'LGALS1', 'S100A16', 'CLIC1', 'MGST1', 'RCAN1', 'TAGLN2', 'NPC2', 'SERPING1', 'C8orf4', 'EMP1', 'APOE', 'CTSB', 'C3', 'LGALS3','EMP3', 'SERPINA3', 'ACTN1', 'PRDX6', 'IGFBP7', 'SERPINE1', 'PLP2', 'MGP', 'CLIC4', 'GFPT2', 'GSN', 'NNMT', 'TUBA1C', 'GJA1', 'TNFRSF1A', 'WWTR1','HILPDA', 'ADM', 'DDIT3', 'NDRG1', 'HERPUD1', 'DNAJB9', 'TRIB3', 'ENO2', 'AKAP12', 'SQSTM1', 'MT1A', 'ATF3', 'NAMPT', 'NRN1', 'SLC2A1','BNIP3', 'LGALS3','INSIG2', 'IGFBP3', 'PPP1R15A', 'VIM', 'PLOD2', 'GBE1', 'SLC2A3', 'FTL', 'WARS', 'ERO1A', 'XPOT', 'HSPA5', 'GDF15', 'ANXA2', 'EPAS1', 'LDHA', 'P4HA1', 'SERTAD1', 'PFKP', 'PGK1', 'EGLN3', 'SLC6A6', 'CA9', 'BNIP3L', 'RPL21', 'TRAM1', 'UFM1', 'ASNS', 'GOLT1B', 'ANGPTL4', 'SLC39A14', 'CDKN1A', 'HSPA9'))
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'MES', search = F)
 # AC module score ####
cd_features <- list(c('CST3', 'S100B', 'SLC1A3', 'HAP1', 'HOPX', 'MT3', 'SPARCL1', 'MLC1', 'GFAP', 'FABP7', 'BCAN', 'PON2', 'METTL7B', 'SPARC', 'GATM', 'RAMP1','PMP2', 'AQP4', 'DBI', 'EDNRB', 'PTPRZ1', 'CLU', 'PMP22', 'ATP1A2', 'S100A16', 'HEY1', 'PCDHGC3', 'TTYH1', 'NDRG2', 'PRCP', 'ATP1B2', 'AGT','PLTP', 'GPM6B', 'F3', 'RAB31', 'PLPP3', 'ANXA5', 'TSPAN7')) 
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'AC', search = F)
 # OPC Module Score ####
cd_features <- list(c('BCAN', 'PLP1', 'GPR17', 'FIBIN', 'LHFPL3', 'OLIG1', 'PSAT1', 'SCRG1', 'OMG', 'APOD', 'SIRT2', 'TNR', 'THY1', 'PHYHIPL', 'SOX2OT', 'NKAIN4', 'PLPPR1', 'PTPRZ1', 'VCAN', 'DBI', 'PMP2', 'CNP', 'TNS3', 'LIMA1', 'CA10', 'PCDHGC3', 'CNTN1', 'SCD5', 'P2RX7', 'CADM2', 'TTYH1', 'FGF12', 'TMEM206', 'NEU4', 'FXYD6', 'RNF13', 'RTKN', 'GPM6B', 'LMF1', 'ALCAM', 'PGRMC1', 'HRASLS', 'BCAS1', 'RAB31', 'PLLP', 'FABP5', 'NLGN3', 'SERINC5', 'EPB41L2', 'GPR37L1'))
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'OPC', search = F)
 # Cell cycle Module Score ####
cd_features <- list(c("UBE2T", "HMGB2", "TYMS", "MAD2L1", "CDK1", "UBE2C", "RRM2", "PBK", "ZWINT","NUSAP1", "PCNA", "BIRC5", "H2AFZ", "FAM64A", "TOP2A", "KIAA0101", "PTTG1","GMNN", "KPNA2", "TUBA1B", "NUF2", "TPX2", "CENPU", "HIST1H4C", "KIF22", "TMPO", "CKS2", "CDCA5", "CENPM", "PRC1", "MCM7", "TMSB15A", "CENPF", "RNASEH2A","RACGAP1", "DUT", "CKS1B", "AURKB", "CCNB2", "DTL", "FEN1", "FANCI", "KIF11","RRM1", "MCM2", "CDC20", "HMGN2", "CCNA2", "TK1", "PKMYT1"))
Adult <- AddModuleScore(object = Adult, features = cd_features, name = 'Cell.Cycle', search = F)

# Mouse
 # NPC1 Module Score ####
cd_features <- list(c("Dll3",'Dll1','Sox4','Tubb3','Hes6','Tagln3','Neu4','Marcksl1','Cd24a','Stmn1','Tcf12','Bex1','Olig1','Map2','Fxyd6','Ptprs','Mllt11','Nppa','Bcan','Mest','Ascl1','Btg2','Dcx','Nxph1','Jpt1','Pfn2','Scg3','Myt1','Chd7','Adgrg1','Tuba1a','Pcbp4','Etv1','Shd','Tnr','Amotl2','Dbn1','Hip1','Abat','Elavl4','Lmf1','Grik2','Serinc5','Tspan13','Elmo1','Glcci1','Sez6l','Lrrn1','Sez6','Sox11'))
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'NPC', search = F)
 # MES1 module score ####
cd_features <- list(c('Chil1', 'Anxa2', 'Anxa1', 'Cd44', 'Vim', 'Mt2', 'C1s1', 'Nampt', 'Efemp1', 'C1ra', 'Sod2', 'Ifitm3', 'Timp1', 'Spp1', 'A2m', 'S100a11','Mt1', 'S100a10', 'Fn1', 'Lgals1', 'S100a16', 'Clic1', 'Mgst1', 'Rcan1', 'Tagln2', 'Npc2', 'Serping1', 'Tcim', 'Emp1', 'Apoe', 'Ctsb', 'C3', 'Lgals3','Emp3', 'Serpina3n', 'Actn1', 'Prdx6', 'Igfbp7', 'Serpine1', 'Plp2', 'Mgp', 'Clic4', 'Gfpt2', 'Gsn', 'Nnmt', 'Tuba1c', 'Gja1', 'Tnfrsf1a', 'Wwtr1'))
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'MES1', search = F)
 # MES2 module score ####
cd_features <- list(c('Hilpda', 'Adm', 'Ddit3', 'Ndrg1', 'Herpud1', 'Dnajb9', 'Trib3', 'Eno2', 'Akap12', 'Sqstm1', 'Mt1', 'Atf3', 'Nampt','Nrn1', 'Slc2a1','Bnip3', 'Lgals3', 'Insig2', 'Igfbp3', 'Ppp1r15a', 'Vim', 'Plod2', 'Gbe1', 'Slc2a3', 'Ftl1', 'Wars','Ero1l', 'Xpot', 'Hspa5', 'Gdf15', 'Anxa2', 'Epas1', 'Ldha', 'P4ha1', 'Sertad1', 'Pfkp', 'Pgk1', 'Egln3', 'Slc6a6', 'Car9', 'Bnip3l', 'Rpl21', 'Tram1', 'Ufm1', 'Asns', 'Golt1b', 'Angptl4', 'Slc39a14', 'Cdkn1a', 'Hspa9'))
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'MES2', search = F)
 # MES1+2 module score ####
cd_features <- list(c('Chil1', 'Anxa2', 'Anxa1', 'Cd44', 'Vim', 'Mt2', 'C1s1', 'Nampt', 'Efemp1', 'C1ra', 'Sod2', 'Ifitm3', 'Timp1', 'Spp1','A2m', 'S100a11','Mt1', 'S100a10', 'Fn1', 'Lgals1', 'S100a16', 'Clic1', 'Mgst1', 'Rcan1', 'Tagln2', 'Npc2', 'Serping1', 'Tcim', 'Emp1', 'Apoe', 'Ctsb', 'C3', 'Lgals3','Emp3', 'Serpina3n', 'Actn1', 'Prdx6', 'Igfbp7', 'Serpine1', 'Plp2', 'Mgp', 'Clic4', 'Gfpt2', 'Gsn', 'Nnmt', 'Tuba1c', 'Gja1', 'Tnfrsf1a', 'Wwtr1','Hilpda', 'Adm', 'Ddit3', 'Ndrg1', 'Herpud1', 'Dnajb9', 'Trib3', 'Eno2', 'Akap12', 'Sqstm1', 'Mt1', 'Atf3', 'Nampt', 'Nrn1', 'Slc2a1','Bnip3', 'Lgals3','Insig2', 'Igfbp3', 'Ppp1r15a', 'Vim', 'Plod2', 'Gbe1', 'Slc2a3', 'Ftl1', 'Wars', 'Ero1l', 'Xpot', 'Hspa5', 'Gdf15', 'Anxa2', 'Epas1', 'Ldha', 'P4ha1', 'Sertad1', 'Pfkp', 'Pgk1', 'Egln3', 'Slc6a6', 'Car9', 'Bnip3l', 'Rpl21', 'Tram1', 'Ufm1', 'Asns', 'Golt1b', 'Angptl4', 'Slc39a14', 'Cdkn1a', 'Hspa9'))
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'MES', search = F)
 # AC module score ####
cd_features <- list(c('Cst3', 'S100b', 'Slc1a3', 'Hap1', 'Hopx', 'Mt3', 'Sparcl1', 'Mlc1', 'Gfap', 'Fabp7', 'Bcan', 'Pon2', 'Mettl7b', 'Sparc', 'Gatm', 'Ramp1','Pmp2', 'Aqp4', 'Dbi', 'Ednrb', 'Ptprz1', 'Clu', 'Pmp22', 'Atp1a2', 'S100a16', 'Hey1', 'Pcdhgc3', 'Ttyh1', 'Ndrg2', 'Prcp', 'Atp1b2', 'Agt','Pltp', 'Gpm6b', 'F3', 'Rab31', 'Plpp3', 'Anxa5', 'Tspan7')) 
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'AC', search = F)
 # OPC Module Score ####
cd_features <- list(c('Bcan', 'Plp1', 'Gpr17', 'Fibin', 'Lhfpl3', 'Olig1', 'Psat1', 'Scrg1', 'Omg', 'Apod', 'Sirt2', 'Tnr', 'Thy1', 'Phyhipl', 'Sox2-ot', 'Nkain4', 'Lppr1', 'Ptprz1', 'Vcan', 'Dbi', 'Pmp2', 'Cnp', 'Tns3', 'Lima1', 'Ca10', 'Pcdhgc3', 'Cntn1', 'Scd5', 'P2rx7', 'Cadm2', 'Ttyh1', 'Fgf12', 'Tmem206', 'Neu4', 'Fxyd6', 'Rnf13', 'Rtkn', 'Gpm6b', 'Lmf1', 'Alcam', 'Pgrmc1', 'Hrasls', 'Bcas1', 'Rab31', 'Pllp', 'Fabp5', 'Nlgn3', 'Serinc5', 'Epb41l2', 'Gpr37l1'))
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'OPC', search = F)
 # Cell cycle Module Score ####
cd_features <- list(c("Ube2t", "Hmgb2", "Tyms", "Mad2l1", "Cdk1", "Ube2c", "Rrm2", "Pbk", "Zwint","Nusap1", "Pcna", "Birc5", "H2afz", "Pimreg", "Top2a", "Pclaf", "Pttg1","Gmnn", "Kpna2", "Tuba1b", "Nuf2", "Tpx2", "Cenpu", "Hist1h4c", "Kif22", "Tmpo", "Cks2", "Cdca5", "Cenpm", "Prc1", "Mcm7", "Tmsb15a", "Cenpf", "Rnaseh2a","Racgap1", "Dut", "Cks1b", "Aurkb", "Ccnb2", "Dtl", "Fen1", "Fanci", "Kif11","Rrm1", "Mcm2", "Cdc20", "Hmgn2", "Ccna2", "Tk1", "Pkmyt1"))
MADR <- AddModuleScore(object = MADR, features = cd_features, name = 'Cell.Cycle', search = F)

##-##########################-##
#### Recover Modules values ####
##-##########################-##
# With % of cells in each region
corners <- rbind(c(1, 1), c(1, -1), c(-1, -1), c(-1, 1)) # 4 variables therefore square
Human <- Adult@meta.data[, c('subtype','NPC1','MES1','AC1','OPC1')]
colnames(Human) <-c('subtype','NPC','MES','AC','OPC')
Human[6:7] <- as.data.frame(as.matrix(Human[2:5])%*%corners)

Human$V1 <- scales::rescale(Human$V1, to=c(-1,1))
Human$V2 <- scales::rescale(Human$V2, to=c(-1,1))

Human <- transform(Human, Phenotype= ifelse(V1<0 & V2>0, 'OPC',
                                            ifelse(V1>0 & V2>0, 'NPC',
                                                   ifelse(V1>0 & V2<0, 'MES', 'AC'))))
Human$subtype <-droplevels(Human$subtype)
Human2 <- as.data.frame.matrix(table(Human[,c(1,8)]))
Human2 <- (Human2/rowSums(Human2))
Human2$subtype<- as.factor(rownames(Human2))
rownames(Human2) <- NULL
Human2$subtype <- factor(Human2$subtype,levels = c("EGFRsnv","EGFRamp, EGFRsnv, PTEN & CDKN2A","EGFRamp + NF1"), ordered = T)
Human2 <- Human2[,c("subtype","AC","MES","NPC","OPC")]
p1 <- ggradar(Human2,  grid.min = 0,
              grid.mid = 0.4,
              grid.max = 0.8,
              fill = T, fill.alpha = 0.2,
              plot.title = 'Human EGFR', legend.title = 'Genomic alterations', 
              group.point.size	=1, group.line.width =0.3,
              background.circle.colour= 'azure',grid.line.width = 0.3,
              gridline.min.colour = "lightgrey", gridline.mid.colour = "royalblue1",gridline.max.colour = "coral",
              label.gridline.min = F,label.gridline.mid = F,label.gridline.max = F)

#########
Mouse <- MADR@meta.data[, c('sample_label','NPC1','MES1','AC1','OPC1')]
Mouse$sample_label <- sub("Egfr glioma 1 hit", "EGFRvIII", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 3 hit female 1", "EGFRvIII + Pten-Cdkn2a-Cas9", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 3 hit male", "EGFRvIII + Pten-Cdkn2a-Cas9", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 4 hit Nf1 tissue", "EGFRvIII + Pten-Cdkn2a-Nf1-Cas9", Mouse$sample_label)
colnames(Mouse) <-c('subtype','NPC','MES','AC','OPC')
corners <- rbind(c(1, 1), c(1, -1), c(-1, -1), c(-1, 1)) # square
Mouse[6:7] <- as.data.frame(as.matrix(Mouse[2:5])%*%corners)

Mouse$V1 <- scales::rescale(Mouse$V1, to=c(-0.8,0.3))
Mouse$V2 <- scales::rescale(Mouse$V2, to=c(-1,0.3))
Mouse <- transform(Mouse, Phenotype= ifelse(V1<0 & V2>0, 'OPC',
                                            ifelse(V1>0 & V2>0, 'NPC',
                                                   ifelse(V1>0 & V2<0, 'MES', 'AC'))))

Mouse$subtype <-droplevels(Mouse$subtype)
Mouse2 <- as.data.frame.matrix(table(Mouse[,c(1,8)]))
Mouse2 <- (Mouse2/rowSums(Mouse2))
Mouse2$subtype<- rownames(Mouse2)
Mouse2 <- Mouse2[c("EGFRvIII","EGFRvIII + Pten-Cdkn2a-Cas9","EGFRvIII + Pten-Cdkn2a-Nf1-Cas9"),c("subtype","AC","MES","NPC","OPC")]
p2 <- ggradar(Mouse2,   values.radar = c("0%", "35%", "70%"),
              grid.min = -0,
              grid.mid = 0.4,
              grid.max = 0.8,
              fill = T, fill.alpha = 0.2, 
              plot.title = 'Mouse EGFR', legend.title = 'Genomic alterations', 
              group.point.size	=1, group.line.width =0.3,
              background.circle.colour= 'azure', grid.line.width = 0.3,
              gridline.min.colour = "lightgrey", gridline.mid.colour = "royalblue1",gridline.max.colour = "coral",
              label.gridline.min = F,label.gridline.mid = F,label.gridline.max = F)

p1/p2

#With actual module values
Human <- Adult@meta.data[, c('subtype','NPC1','MES1','AC1','OPC1')]
colnames(Human) <-c('subtype','NPC','MES','AC','OPC')
Human$subtype <-droplevels(Human$subtype)

MaxValue <-as.data.frame(c('Max', (Human %>% summarise(across(colnames(Human[2:5]), max),.groups = 'drop')  %>%  as.data.frame())))
colnames(MaxValue) <- colnames(Human)
percent99 <- as.data.frame(Human %>% do(data.frame(subtype = '99 percentile',
                                                   NPC=quantile(.$NPC, probs=c(0.99)),
                                                   MES=quantile(.$MES, probs=c(0.99)),
                                                   AC=quantile(.$AC, probs=c(0.99)),
                                                   OPC=quantile(.$OPC, probs=c(0.99)))))
percent01 <- as.data.frame(Human %>% do(data.frame(subtype = '01 percentile',
                                                   NPC=quantile(.$NPC, probs=c(0.01)),
                                                   MES=quantile(.$MES, probs=c(0.01)),
                                                   AC=quantile(.$AC, probs=c(0.01)),
                                                   OPC=quantile(.$OPC, probs=c(0.01)))))
MinValue <-as.data.frame(c('Min', (Human %>% summarise(across(colnames(Human[2:5]), min),.groups = 'drop')  %>%  as.data.frame())))
colnames(MinValue) <- colnames(Human)

Human2 <- Human %>% group_by(subtype) %>% 
  summarise(across(colnames(Human[2:5]), mean),.groups = 'drop')  %>%
  as.data.frame()
Human2 <- rbind(Human2, percent99,percent01)
Human2 <- rbind(Human2,MaxValue,MinValue)
Human2[,2:5] <- scale(Human2[,2:5])
#Human2 <- scale(Human2)
Human2 <- Human2[1:3,]
Human2$subtype <- factor(Human2$subtype,levels = c("EGFRsnv","EGFRamp, EGFRsnv, PTEN & CDKN2A","EGFRamp + NF1"), ordered = T)
Human2 <- Human2[,c("subtype","AC","MES","NPC","OPC")]

#########
Mouse <- MADR@meta.data[, c('sample_label','NPC1','MES1','AC1','OPC1')]
Mouse$sample_label <- sub("Egfr glioma 1 hit", "EGFRvIII", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 3 hit female 1", "EGFRvIII + Pten-Cdkn2a-Cas9", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 3 hit male", "EGFRvIII + Pten-Cdkn2a-Cas9", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 4 hit Nf1 tissue", "EGFRvIII + Pten-Cdkn2a-Nf1-Cas9", Mouse$sample_label)

colnames(Mouse) <-c('subtype','NPC','MES','AC','OPC')
Mouse$subtype <-as.factor(Mouse$subtype)

MaxValue <-as.data.frame(c('Max', (Mouse %>% summarise(across(colnames(Mouse[2:5]), max),.groups = 'drop')  %>%  as.data.frame())))
colnames(MaxValue) <- colnames(Mouse)
percent99 <- as.data.frame(Mouse %>% do(data.frame(subtype = '99 percentile',NPC=quantile(.$NPC, probs=c(0.99)),
                                                   MES=quantile(.$MES, probs=c(0.99)),AC=quantile(.$AC, probs=c(0.99)),
                                                   OPC=quantile(.$OPC, probs=c(0.99)))))
percent01 <- as.data.frame(Mouse %>% do(data.frame(subtype = '01 percentile',NPC=quantile(.$NPC, probs=c(0.01)),
                                                   MES=quantile(.$MES, probs=c(0.01)),AC=quantile(.$AC, probs=c(0.01)),
                                                   OPC=quantile(.$OPC, probs=c(0.01)))))
MinValue <-as.data.frame(c('Min', (Mouse %>% summarise(across(colnames(Mouse[2:5]), min),.groups = 'drop')  %>%  as.data.frame())))
colnames(MinValue) <- colnames(Mouse)
Mouse2 <- Mouse %>% group_by(subtype) %>% 
  summarise(across(colnames(Mouse[2:5]), mean),.groups = 'drop')  %>%
  as.data.frame()
Mouse2 <- rbind(Mouse2, percent99,percent01)
Mouse2 <- rbind(Mouse2, MaxValue,MinValue)
Mouse2[,2:5] <- scale(Mouse2[,2:5])
Mouse2 <- Mouse2[1:3,]
Mouse2 <- Mouse2[,c("subtype","AC","MES","NPC","OPC")]
##### plots

p1 <- ggradar(Human2,   grid.min = min(Human2[,2:5])-.1,
              grid.mid = mean(c(min(Human2[,2:5]), max(Human2[,2:5]))),
              grid.max = max(Human2[,2:5])+.05,
              centre.y =  min(Human2[,2:5])-0.25,
              fill = T, fill.alpha = 0.1,
              plot.title = 'Human EGFR', legend.title = 'Genomic alterations', 
              group.point.size	=1.5, group.line.width =0.5,
              background.circle.colour= 'azure',grid.line.width = 0.3,
              gridline.min.colour = "lightgrey", gridline.mid.colour = "#B2BEB5",gridline.max.colour = "#36454F",
              label.gridline.min = F,label.gridline.mid = F,label.gridline.max = F,
              group.colours = c('#1fe98b','#5ba1d3','#fd9f4b'))+guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.key.width= unit(0.6, 'cm'))



p2 <- ggradar(Mouse2,  
              grid.min = min(Mouse2[,2:5])-.1,
              grid.mid = mean(c(min(Mouse2[,2:5]), max(Mouse2[,2:5]))),
              grid.max = max(Mouse2[,2:5])+.05,
              centre.y =  min(Mouse2[,2:5])-0.25,
              fill = T, fill.alpha = 0.1,
              plot.title = 'Mouse EGFR', legend.title = 'Genomic alterations', 
              group.point.size	=1.5, group.line.width =0.5,
              background.circle.colour= 'azure',grid.line.width = 0.3,
              gridline.min.colour = "lightgrey", gridline.mid.colour = "#B2BEB5",gridline.max.colour = "#36454F",
              label.gridline.min = F,label.gridline.mid = F,label.gridline.max = F,
              group.colours=c('#2ba02b','#1f77b4','#ff7f0f'))+guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.key.width= unit(0.6, 'cm'))

p1/p2

# Combining non-Nf1 subtypes #####
Human <- Adult@meta.data[, c('subtype','NPC1','MES1','AC1','OPC1')]
colnames(Human) <-c('subtype','NPC','MES','AC','OPC')
Human$subtype <-droplevels(Human$subtype)
Human$subtype2 <- with(Human, ifelse(subtype == 'EGFRamp + NF1', 'NF1', 'Non-NF1'))

#MaxValue <-as.data.frame(c('Max', (Human %>% summarise(across(colnames(Human[2:5]), max),.groups = 'drop')  %>%  as.data.frame())))
#colnames(MaxValue) <- colnames(Human[c(6,2,3,4,5)])
percent99 <- as.data.frame(Human %>% do(data.frame(subtype2 = '99 percentile',
                                                   NPC=quantile(.$NPC, probs=c(0.99)),
                                                   MES=quantile(.$MES, probs=c(0.99)),
                                                   AC=quantile(.$AC, probs=c(0.99)),
                                                   OPC=quantile(.$OPC, probs=c(0.99)))))
percent01 <- as.data.frame(Human %>% do(data.frame(subtype2 = '01 percentile',
                                                   NPC=quantile(.$NPC, probs=c(0.01)),
                                                   MES=quantile(.$MES, probs=c(0.01)),
                                                   AC=quantile(.$AC, probs=c(0.01)),
                                                   OPC=quantile(.$OPC, probs=c(0.01)))))
#MinValue <-as.data.frame(c('Min', (Human %>% summarise(across(colnames(Human[2:5]), min),.groups = 'drop')  %>%  as.data.frame())))
#colnames(MinValue) <- colnames(Human[c(6,2,3,4,5)])

Human2 <- Human %>% group_by(subtype2) %>% 
  summarise(across(colnames(Human[2:5]), mean),.groups = 'drop')  %>%
  as.data.frame()
Human2 <- rbind(Human2, percent99,percent01)
#Human2 <- rbind(Human2,MaxValue,MinValue)
Human2[,2:5] <- scale(Human2[,2:5])
Human2 <- Human2[1:2,]
Human2 <- Human2[,c("subtype2","AC","MES","NPC","OPC")]

#########
Mouse <- MADR@meta.data[, c('sample_label','NPC1','MES1','AC1','OPC1')]
Mouse$sample_label <- sub("Egfr glioma 1 hit", "EGFRvIII", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 3 hit female 1", "EGFRvIII + Pten-Cdkn2a-Cas9", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 3 hit male", "EGFRvIII + Pten-Cdkn2a-Cas9", Mouse$sample_label)
Mouse$sample_label <- sub("Egfr glioma 4 hit Nf1 tissue", "EGFRvIII + Pten-Cdkn2a-Nf1-Cas9", Mouse$sample_label)

colnames(Mouse) <-c('subtype','NPC','MES','AC','OPC')
Mouse$subtype <-as.factor(Mouse$subtype)
Mouse$subtype2 <- with(Mouse, ifelse(subtype == 'EGFRvIII + Pten-Cdkn2a-Nf1-Cas9', 'Nf1', 'Non-Nf1'))

#MaxValue <-as.data.frame(c('Max', (Mouse %>% summarise(across(colnames(Mouse[2:5]), max),.groups = 'drop')  %>%  as.data.frame())))
#colnames(MaxValue) <- colnames(Mouse[c(6,2,3,4,5)])
percent99 <- as.data.frame(Mouse %>% do(data.frame(subtype2 = '99 percentile',NPC=quantile(.$NPC, probs=c(0.99)),
                                                   MES=quantile(.$MES, probs=c(0.99)),AC=quantile(.$AC, probs=c(0.99)),
                                                   OPC=quantile(.$OPC, probs=c(0.99)))))
percent01 <- as.data.frame(Mouse %>% do(data.frame(subtype2 = '01 percentile',NPC=quantile(.$NPC, probs=c(0.01)),
                                                   MES=quantile(.$MES, probs=c(0.01)),AC=quantile(.$AC, probs=c(0.01)),
                                                   OPC=quantile(.$OPC, probs=c(0.01)))))
#MinValue <-as.data.frame(c('Min', (Mouse %>% summarise(across(colnames(Mouse[2:5]), min),.groups = 'drop')  %>%  as.data.frame())))
#colnames(MinValue) <- colnames(Mouse[c(6,2,3,4,5)])
Mouse2 <- Mouse %>% group_by(subtype2) %>% 
  summarise(across(colnames(Mouse[2:5]), mean),.groups = 'drop')  %>%
  as.data.frame()
Mouse2 <- rbind(Mouse2, percent99,percent01)
Mouse2 <- rbind(Mouse2, MaxValue,MinValue)
Mouse2[,2:5] <- scale(Mouse2[,2:5])
Mouse2 <- Mouse2[1:2,]
Mouse2 <- Mouse2[,c("subtype2","AC","MES","NPC","OPC")]
##### plots

p1<- ggradar(Human2,
             grid.min = min(Human2[,2:5])-(0.2*(max(Human2[,2:5])- min(Human2[,2:5]))),
             grid.mid = mean(c(min(Human2[,2:5]), max(Human2[,2:5]))),
             grid.max = max(Human2[,2:5])+(0.2*(max(Human2[,2:5])- min(Human2[,2:5]))),
             centre.y =  min(Human2[,2:5])-(0.4*(max(Human2[,2:5])- min(Human2[,2:5]))),
             fill = T, fill.alpha = 0.1,
             plot.title = 'Human EGFR', legend.title = 'Genomic alterations', 
             group.point.size	=1.5, group.line.width =0.5,
             background.circle.colour= 'white',grid.line.width = 0.3,
             gridline.min.colour = "lightgrey", gridline.mid.colour = "#B2BEB5",gridline.max.colour = "#36454F",
             label.gridline.min = F,label.gridline.mid = F,label.gridline.max = F,
             group.colours = c('#fd9f4b','#1f77b4'))+guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.key.width= unit(0.6, 'cm'))



p2 <- ggradar(Mouse2,  
              grid.min = min(Mouse2[,2:5])-(0.2*(max(Mouse2[,2:5])- min(Mouse2[,2:5]))),
              grid.mid = mean(c(min(Mouse2[,2:5]), max(Mouse2[,2:5]))),
              grid.max = max(Mouse2[,2:5])+(0.2*(max(Mouse2[,2:5])- min(Mouse2[,2:5]))),
              centre.y =  min(Mouse2[,2:5])-(0.4*(max(Mouse2[,2:5])- min(Mouse2[,2:5]))),
              fill = T, fill.alpha = 0.1,
              plot.title = 'Mouse EGFR', legend.title = 'Genomic alterations', 
              group.point.size	=1.5, group.line.width =0.5,
              background.circle.colour= 'white',grid.line.width = 0.3,
              gridline.min.colour = "lightgrey", gridline.mid.colour = "#B2BEB5",gridline.max.colour = "#36454F",
              label.gridline.min = F,label.gridline.mid = F,label.gridline.max = F,
              group.colours=c('#fd9f4b','#1f77b4'))+guides(color = guide_legend(override.aes = list(size = 2)))+
  theme(legend.key.width= unit(0.6, 'cm'))

p1/p2




