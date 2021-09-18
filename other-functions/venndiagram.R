
elasticnet_genes = c("PGK1", "IMP5", "SHCBP1", "FIBCD1", "LOC100128977", "PPFIA3", "WNT3A", "PCYT1A", "CRISP3", "TPT1", "DCTPP1", "IRF2", "PSD2", "DSG1", "FAM9C", "PELO", "XRCC4", "LSG1", "SERPINA3", "GCET2", "IL18", "IYD", "LOC220729", "SLC20A2", "PCMT1", "TCP1", "HSPA8", "KCNJ13", "RPL3", "LHFPL3", "JAK1", "ABCG4", "MED17", "SFRP5", "GLUL", "BAMBI", "CSN3", "C17orf64", "NIPA2", "ZNF485", "NDST4", "MURC", "DHX16", "AKR1E2", "C8orf55", "OPN4", "SPINT1", "C11orf20", "GPR172A", "XG")
hubcox_genes = c("LOC100128977", "SHCBP1", "PPFIA3", "PGK1", "PCYT1A", "IMP5", "TPT1", "IRF2", "FIBCD1", " DCTPP1", "GLUL", "TCP1", "AKR1E2", "MED17", "GPR172A", "LSG1", "LOC220729", "JAK1", "HSPA8", "PELO", "XRCC4", "RPL3", "C8orf55", "PCMT1", "GCET2", "IL18", "PTGES3", "SERPINA3", "NIPA2", "OPA1", "DHX16", "PWP2", "XG", "IYD", "ZNF485", "SPINT1", "CLDN7", "HSPA9", "WNT3A", "PSME1", "TANK", "RUNX1", "ABCG4", "GFI1")
orphancox_genes = c("CRISP3", "PPFIA3", "WNT3A", "PGK1", "IMP5", "LOC100128977", "FIBCD1", "SHCBP1", "PCYT1A", "PSD2", "SERPINA3", "IRF2", "LHFPL3", "FAM9C", "C17orf64", "GCET2", "LSG1", "TPT1", "SLC20A2", "PELO", "IYD", "PCMT1", "BAMBI", "DCTPP1", "SFRP5", "XRCC4", "KCNJ13", "CSN3", "DSG1", "RPL3", "IL18", "LOC220729", "NDST4", "ABCG4", "AKR1E2", "MURC", "CRISP2", "MED17", "JAK1", "MAFA", "GLUL", "HSPA8", "TCP1", "C11orf20", "DHX16")
tcox_brca_genes = c("VHL", "HYI", "MIER1", "MRE11A", "ACAP2", "MCTS1", "PSME2", "ZNF32", "PCGF5", "RPL14", "ARID1B", "DCTPP1", "SLC35C2", "RPL29", "APOOL", "COMTD1", "UBXN7", "SECISBP2L", "JAK1", "GPR172A", "GPR107", "PSENEN", "ZFC3H1", "PCYT1A", "PPFIA3", "CHUK", "DIP2B", "FAM98B", "SHARPIN", "SDR39U1", "BOLA1", "DGAT1", "COMMD5", "STXBP5", "RSPRY1", "GNB1L", "ZFYVE19", "PSME1", "RPL3", "LOC728323", "NQO2", "SPINT1", "C11orf20", "PFKL", "PTPMT1", "C21orf33", "TMEM60")
tcox_prad_genes = c("C9orf80", "QTRTD1", "HN1", "SLC23A3", "REX04", "POLE3", "C14orf118", "OXSR1", "ITGA2B", "C19orf76", "SURF6", "KIAA1958", "C9orf102", "SYNJ1", "POLD1", "COMT", "ZFP37", "LARP7", "TNFRSF8", "ZNF749", "AKAP8L", "TIAL1", "YY1", "NAPSA", "HELLS", "TSEN54", "C18orf25", "RPE", "ZNF189", "THRSP", "WDR82", "FAM169B", "LHX5", "ZNF165", "AGAP6", "RANBP6", "LOC100272146", "BBS7", "KBTBD8", "C14orf129", "ESCO1", "ZNF335", "MRPL30", "TEKT5", "ZNF330", "ZNF70", "GTF2H1")
rsf_genes = c("CEACAM5", "STXBP5", "ENC1", "DNAJC22", "GLT25D2", "CEL", "SCG5", "HPDL", "CALML5", "XRCC4", "FIBCD1", "WNT3A", "PHB", "HSPA8", "DUS1L", "RBBP8", "CCDC54", "C6orf141", "GPR37L1", "ZNF654", "DAB2", "CEACAM1", "MLLT6", "ADAMTSL1", "ADAMTS7", "TNFSF4", "LOC728323", "DNAJC14", "NACC1", "ABHD10", "GADL1", "NLE1", "MYO16", "OPLAH", "ZNF707", "SPRY4", "MGC21881", "CA11", "MTMR7", "NXN", "MRPL38", "BCLAF1", "ZCCHC9", "CHERP", "PLXNB2", "PUF60", "TSGA10", "LHX5", "BAZ2A", "C21orf57")
tcox_all_genes = c("VHL", "HYI", "MIER1", "MRE11A", "ACAP2", "MCTS1", "PSME2", "ZNF32", "PCGF5", "RPL14", "ARID1B", "DCTPP1", "SLC35C2", "RPL29", "APOOL", "COMTD1", "UBXN7", "SECISBP2L", "JAK1", "GPR172A", "GPR107", "PSENEN", "ZFC3H1", "PCYT1A", "PPFIA3", "CHUK", "DIP2B", "FAM98B", "SHARPIN", "SDR39U1", "BOLA1", "DGAT1", "COMMD5", "STXBP5", "RSPRY1", "GNB1L", "ZFYVE19", "PSME1", "RPL3", "LOC728323", "NQO2", "SPINT1", "C11orf20", "PFKL", "PTPMT1", "C21orf33", "TMEM60", "C9orf80", "QTRTD1", "HN1", "SLC23A3", "REX04", "POLE3", "C14orf118", "OXSR1", "ITGA2B", "C19orf76", "SURF6", "KIAA1958", "C9orf102", "SYNJ1", "POLD1", "COMT", "ZFP37", "LARP7", "TNFRSF8", "ZNF749", "AKAP8L", "TIAL1", "YY1", "NAPSA", "HELLS", "TSEN54", "C18orf25", "RPE", "ZNF189", "THRSP", "WDR82", "FAM169B", "LHX5", "ZNF165", "AGAP6", "RANBP6", "LOC100272146", "BBS7", "KBTBD8", "C14orf129", "ESCO1", "ZNF335", "MRPL30", "TEKT5", "ZNF330", "ZNF70", "GTF2H1")

# Venn diagram
library("VennDiagram")

setwd("~/Documents/msc-beatriz-r-correia/")

# venn_diagram1 <- venn.diagram( x = list("TCox-BRCA" = tcox_brca_genes,
#                                         "TCox-PRAD" = tcox_prad_genes),
#                               imagetype="png", filename = NULL, cat.fontface="bold", cat.cex=1.4, palette="ggplot2",
#                               label.col = c("black", "black", "black"),
#                               fontface="bold",cex = 1.8, fill = c("#bfbd95", "#faf593"))
# grid.draw(venn_diagram1)

venn_diagram2 <- venn.diagram( x = list( "1" = elasticnet_genes,
                                        "2" = hubcox_genes,
                                        "3" = orphancox_genes,
                                        "4" = tcox_all_genes,
                                        "5" = rsf_genes),
                              imagetype="png", filename = NULL, cat.fontface="bold", cat.cex=1.4, palette="ggplot2",
                              label.col = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
                              fontface="bold",cex = 1.8, fill = c("#96ebcd", "#a8b0f0", "#ff9c6b", "#fff673", "#ffabff"))
grid.draw(venn_diagram2)

# venn_diagram3 <- venn.diagram( x = list( "1" = elasticnet_genes,
#                                          "2" = hubcox_genes,
#                                          "3" = orphancox_genes,
#                                          "4" = tcox_all_genes),
#                                imagetype="png", filename = NULL, cat.fontface="bold", cat.cex=1.4, palette="ggplot2",
#                                label.col = c("black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black", "black"),
#                                fontface="bold",cex = 1.8, fill = c("#96ebcd", "#a8b0f0", "#ff9c6b", "#fff673"))
# grid.draw(venn_diagram3)