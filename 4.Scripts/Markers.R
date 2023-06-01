## Markers

#Markers------------------------------------------------------------------------
#Taken from Fiorenzano et al. Nature Comm 2021, Sozzi et al.2022, 
#Kanton et al. Nature 2019, Birtele et al.

Pluripotency.M <- c("NANOG","POU5F1","BMP4","WNT3A","WNT5A","LEF1","TCF7","AXIN2","DPPA3","OTX1","TBX3","AFP","CLDN1")
Mesoderm.M <- c("TBXT","NODAL","EOMES","GATA3","GATA4","PECAM1","HAND1","MSX2","CER1")
Endoderm.M <- c("SOX17","FOXA2","GATA6","SERPINA1","VAV3","PAX3","ARID5B","MIXL1","CDX2")
Neuroectoderm.M <- c("FOXG1","TUBB3","NCAM1","MAP2","NOTCH1","RBFOX3","SIX3","VAX1","CACNG8","TENM1","UNC5A","FABP7","PPP2R2B","DLK1","PCDHAC2","GLDN","GBX2","NECAB2","LAMP5","LINGO1","JPH4","SLC8A2","EPHA2")

NeuralStemCells.M <- c("TOP2A","MKI67","CDK1","CLSPN","CENPF","PTTG1","ASCL1","PAX6","SOX2")
RadialGlia.M <- c("VIM","PLP1","EDNRB","GFAP","FABP7","MKI67","HMGB2","ASCL1","TFDP2","PAX6","SOX2","HES1","EOMES")
OuterRadialGlia.M <-c("HOPX","PTPRZ1","TNC","FAM107A","MOXD1")
Neurogenesis.M <- c("SOX2","CENPF","HES1","HES5","NES","FABP7","ASCL1")
FloorPlateProgenitors.M <- c("FABP7","HES5","MSX1","SHH","CORIN","FOXA2","LMX1B","OTX2","SOX6","SLC1A3","SOX2","WNT5A")
Neuroblast.M <- c("NEUROD2","NCAM1","SYT1","SLC17A6","SLC32A1","GAD2","PAX6","HES1","EOMES")
NeuronPrecursors.M <- c("SYT1","STMN2","DCX","NCAM1")
Neurons.M <- c("GAP43","MAP2","NNAT","DCX","STMN2","GRIA2","DLX5","NCAM1","NEUROD6","NEUROD2")
MatureNeurons.M <- c("SYT1","NSG2","SNAP25","DLG4","SYN1","NEUN","MAP2","SYP","GRIN2A","GRIA","GRIN2B","GAP43")
VentralForebrainProgenitors.M <- c("PAX6","DLX1","DLX2","DLX5")
HindbrainProgenitors.M <- c("LHX9","GBX2")

OcularTroclearNucleus.M <- c("ISL1","PHOX2B","NEFL","NEFM","LINC00682","SYT12","GAD2")
InhibitoryNeurons.M <- c("GAD1","GAD2","GABRA1","GABRA2","NEFL","SST","ISL1","SLC32A1","DLX1","DLX2","DLX5","THY1","CCK")
ExcitatoryNeurons.M <- c("GRIN2B","BCL11B","SLC1A2","SIX3","NEUROD6","SLC17A7","GNB1","SATB2","NEUROD2","THY1","EMX1","TBR1")
CorticalEN.M <- c("MAP2","ENO2","FOXG1","NEUROD2","NEUROD6","SLC17A7","TBR1","SATB2","BCL11B","RELN","CUX1","BRN2")
MedialandCaudalGanglionicEminenceIN.M <- c("MAP2","ENO2","FOXG1","GAD1","GAD2","SLC32A1","DLX1","DLX2","DLX5")
LateralGanglionicEminenceIN.M <- c("MAP2","ENO2","FOXG1","DLX2","DLX5","GAD2","SLC32A1","ISL1","EBF1")
DiencephalicEN.M <- c("MAP2","ENO2","EBF1","NHLH2","SLC17A6","LHX9","GBX2","SHOX2")
DiencephalicIN.M <- c("MAP2","ENO2","NHLH2","PCP4","RSPO3","RELN")
MesencephalicandRhombencephalicEN.M <- c("MAP2","ENO2","EBF1","NHLH2","SLC17A6","LHX9","TFAP2B")
MesencephalicandRhombencephalicIN.M <- c("MAP2","ENO2","GAD1","GAD2","SLC32A1","GATA3","OTX2","SOX14")
GlutamatergicNeurons.M <- c("SLC17A6","SLC17A7","SLC17A8","SLC1A1","SLC1A2","SLC1A6")
GabaergicNeurons.M <-c("GAD1","GAD2","SLC32A1","SLC6A1")
DopaminergicNeurons.M <-c("TH","DDC","SLC18A2","SLC18A3","SLC6A3")
NoradrenergicNeurons.M <- c("TH","DDC","DBH","SLC6A2","SLC18A2","SLC18A3","SLC39A11","SLC9B2")
AdrenergicNeurons.M <- c("TH","DDC","DBH","PNMT","SLC18A2","NPFF","SLC12A7")
GlicinergicNeurons.M <-c("SLC6A9","SLC32A1")
SerotoninergicNeurons.M <- c("SLC6A4","TPH2","SLC18A2","DDC","FEV","SLC22A3","ESM1","GRIA2","SLC17A8","GATA3","SYT12","TPH1")
CholinergicNeurons.M <-c("SLC5A7","CHAT","ACHE","SLC18A3")

Astrocytes.M <- c("AQP4","GFAP","EDNRB","GJA1","SOX9","SLC1A3","RFX4","S100B","AGT")
OPC.M <- c("OLIG1","OLIG2","NKX2-2","SOX10","PDGFRA","PLP1","MPZ","S100B","MBP","CSPG4")
Oligo.M <- c("OLIG1","OLIG2","S100B","MOG","MBP")
VLMCs.M <- c("COL1A1","COL1A2","DCN","PDGFRA","FBLN1","FBLN2","S100B","S100A11","LUM","IFITM2","PCOLCE","OLFML3","RBP1","CPXM1","MFAP2","MMP2","ISLR2","ADGRD1")
Mural.M <- c("COL1A1","COL1A2","PDGFRA","PDGFRB","IFITM1","FBLN1","S100A11")
RedBloodCells.M <- c("SOD1","BPGM","HMBS","SMIM1")
PericytesEndotelial.M <- c("CLDN5","RGS5","LAMC1","SLC2A1","GNG11","KDR","VIM", "FLT1")
Microglia.M <- c("CD74","CD84","TREM2","ACY3","TYROBP","CLDN5","AIF1","ITGAM")
ChoroidPlexus.M <- c("BMP5","MSX2","CHMP1A","PRLR","SLC29A4","FOXJ1","MSX1","HEY1","SLC26A11","SLC23A2","SLC16A2","SLC4A10")

DAneurogenesis.M <- c("TH","EN1","PBX1","STMN2","SHH","CORIN","SOX6","FOXA2","LMX1A","LMX1B","DDC","NR4A2","KCNJ6","SLC6A3")
DAneurons.M <- c("FOXA2","LMX1A","SLC6A3","TH","NR4A2","PBX1")
DAneuronsHighTH.M <- c("TH","NR4A2","DDC","SLC18A2","SLC6A3","SLITRK5","KCNJ6","ROBO2","LSAMP","NTM","JPH4","SNTG1")
DAneuronsLowTH.M <- c("ASCL1","NES","SOX2","NEUROG2","GNG5","NHLH1","NFIA","NFIB","RND3","CHD7","RFX4","HES6")
DAprog1.M <- c("SOX6","SOX2","ANXA1","NES","LPL","ZDHHC2","ID1","ID2","NKX2-2","CCND1","CENPF","E2F3")
DAprog2.M <- c("NEAT1","SNCG","CDH2","NHLH2")
DAneurons1.M <- c("PBX1","EN1","NR4A2","DDC","FOXA2","SEZ6L2","TH","DLK1","POU2F2","GFRA1","SEZ6L2.1","NONF","C1QL1","CXCR4","OTX2","SLC17A6","MT3","RELN")
DAneurons2.M <- c("PBX1","EN1","NR4A2","DDC","FOXA2","SEZ6L2","TH","DLK1","POU2F2","GFRA1","SEZ6L2.1","CCK","CBLN1","COX17","NEFM","CALB1","TACRS","CD24","GRP","SLC6A3","SLC18A2","ALDH1A1","NETO2","SYT1")
DA.L1.M <- c("CAMK2N1","SLC2A1","SLC1A6","CTXN1","SLC16A3","NKX6-1","EBF2","MEG3","HES6","RBFOX2","NKAIN4","CXCR4","RTN1","DCX","OLFM1","IGDCC3","NNAT","NEFL","NEFM","NHLH2","SNCA","CD24","SNCG","ZDHHC2","DCC","DLK1","SLC6A3","KCNJ6")
DA.L2.M <- c("EN2","MPC1","SRP54","ANXA1","ANXA2","SYAP1","TOMM7","ATP5ME","NDUFA1","COX6A1","COX7A2","COX17")
DA.L3.M <- c("PCP4","GPM6A","STMN2","NRXN2","SLC5A3","SLC8A1","SYT13","CAMK2N2","SEZ6L","SEZ6L2","CACNA2D1","SDC2","FJX1","NRIP3","LMX1A","FOXA2","POU2F2","ID4","CALB2","EN1","GRP","SLC18A2","PLXNC1","NR4A2","LMO3","SOX6","OTX2","CALB1","FOXP1")

SubstantiaNigraDevelopment.GO.M <- c("TTBK1","RHOA","ATP5PF","HSPA5","DYNLL1","CNP","YWHAQ","SEC16A","RAD1","ZNF148","PADI2","CDC42","CALM1","CALM2","MBP","ACTB","FGF2","S100A1","CASP5","FGF9")

ActivityDependentGenes.M <- c("")

#From Poulin 2014 DA diversity -------------------------------------------------
Housekeeping.M <- c("ACTB","GAPDH","HPRT")
PDlinked.M <- c("ATP13A2", "LRRK2", "PARK2", "PARK7", "PINK1", "SNCA")
DAneurons.Poulin.M <- c("DDC","TH","SLC6A3","SLC18A2")
DAtranscriptionfactors <- c("FOXA1","FOXA2","EN1","EN2","LMX1B","PITX3","NR4A2")
SNcenriched.M <- c("SOX6","SLC6A3","ZDHHC2","NRIP3","CHRNA4","KCNJ6","CLSTN2","FAM184A","KCNS3","SATB1","DRD2","RABC3","NDNF","SEMA6D","SNCG","GUCY2C","IGF1","FGF1","PITX3")
VTAenriched.M <- c("FOXA2","EFNB3","CALB2","CCK","SLC17A6","CALB1","LMX1A","VGF","EN2","GRP","FZD1","SDC2","LPL","B4GALT1","EN1","MARCKS","VIP","CACNA1D","OTX2","NEFM","ADCYAP1","TACR3")
DAdiversitygenes.M <- c("SOX6","SNCG","NDNF","IGF1","FOXA2","LMX1A","ALDH1A1","SLC32A1","SATB1","CLSTN2","ADCYAP1","LPL","OTX2","VIP","CHRNA4","GSG1L","SMCA","NTF3")


#Colors-------------------------------------------------------------------------

TwoColors <- c("#7AAABA","#EECF82") 
DotplotColors <- c("white","#245D70") 
UMAPColors7 <- c("#80C3C6","#CF6EAB","#4391BA","#122E7B","#255838","#EB8C44","#60A762")
UMAPColors9 <- c("#EB8C44","#80C3C6","#4391BA","#122E7B","#CF6EAB","#E8312D","#60A762","#255838","#9923A3")

#GO terms-----------------------------------------------------------------------

Glycolysis.GO <- "GO:0006096" # Glycolysis
ERstress.GO <- "GO:0034976" # ER-stress

Gliogenesis.GO <- "GO:0042063" # Gliogenesis
NeuronMaturation.GO <- "GO:0042551" #Neuron maturation
NeuronDifferentiation.GO <- "GO:0030182" #Neuron differentiation
NeuronMigration.GO <- "GO:0001764" #Neuron migration
NeuronDevelopment.GO <- "GO:0048666" #Neuron development
NeuronProjection.GO <- "GO:0043005" #Neuron projection
NeuronDeath.GO <- "GO:0070997" #Neuron death
NeuronFateDetermination.GO <- "GO:0048664" #Neuron Fate Determination

CellularResponseToHypoxia.GO <- "GO:0071456" #Cellular Response To Hypoxia
DetectionOfHypoxia.GO <- "GO:0070483" #Detection of hypoxia
HIF1AsignalingPathway.GO <- "GO:0097411" #Hypoxia-inducible factor-1alpha signaling pathway
ResponseToHypoxia.GO <- "GO:0001666" #Response to hypoxia
IntrinsicApoptoticSignalingPathwayInResponseToHypoxia.GO <- "GO:1990144" #Intrinsic apoptotic signaling pathway in response to hypoxia

ApoptoticProcess.GO <- "GO:0006915" #Apoptotic process
RegulationOfApoptoticProcess.GO <- "GO:0042981" #Regulation of apoptotic process
NegativeRegulationOfApoptoticProcess.go <- "GO:0043066" #Negative regulation of apoptotic process
ApoptoticNuclearChanges.GO <- "GO:0030262" #Apoptotic nuclear changes
ApoptoticChromosomeCondensation.GO <- "GO:0030263" #Apoptotic chromosome condensation
ApoptoticDNAFragmentation.GO <- "GO:0006309" #Apoptotic DNA fragmentation
ApoptoticMitochondrialChanges.GO <- "GO:0008637" #Apoptotic mitochondrial changes
ApoptoticCellClearance.GO <- "GO:0043277" #Apoptotic cell clearance
Apoptosome.GO <- "GO:0043293" #Apoptosome
MitochondrialFragmentationInvolvedInApoptoticProcess.GO <- "GO:0043653" #Mitochondrial fragmentation involved in apoptotic process
RecognitionOfApoptoticCell.GO <- "GO:0043654" #Recognition of apoptotic cell
CellularComponentDisassemblyInvolvedInExecutionPhaseOfApoptosis.GO <-"GO:0006921" #Cellular component disassembly involved in execution phase of apoptosis

IntrinsicApoptoticSignalingPathwayInResponseToOsmoticStress <- "GO:0008627" #Intrinsic apoptotic signaling pathway in response to osmotic stress
IntrinsicApoptoticSignalingPathwayInResponseToHypoxia.GO <- "GO:1990144" #Intrinsic apoptotic signaling pathway in response to hypoxia
IntrinsicApoptoticSignalingPathwayByp53ClassMediator.GO <- "GO:0072332" #Intrinsic apoptotic signaling pathway by p53 class mediator
IntrinsicApoptoticSignalingPathwayInResponseToDNADamageByp53ClassMediator.GO <- "GO:0042771" #Intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator
IntrinsicApoptoticSignalingPathwayInResponseToNitrosativeStress.GO <- "GO:1990442" #Intrinsic apoptotic signaling pathway in response to nitrosative stress
IntrinsicApoptoticSignalingPathwayInResponseToDNADamage.GO <- "GO:0008630" #Intrinsic apoptotic signaling pathway in response to DNA damage
IntrinsicApoptoticSignalingPathwayInResponseToOxidativeStress.GO <- "GO:0008631" #Intrinsic apoptotic signaling pathway in response to oxidative stress

ExtrinsicApoptoticSignalingPathwayViaDeathDomainReceptors.GO <- "GO:0008625" #Extrinsic apoptotic signaling pathway via death domain receptors

GlialCellApoptoticProcess.GO <- "GO:0034349" #Glial cell apoptotic process
RegulationOfNeuronApoptoticProcess.GO <- "GO:0043523" #Regulation of neuron apoptotic process

#Wnt pathway
WntPathway <- "GO:0016055"
RegulationOfWntSignalingPathway <- "GO:0030111"
CanonicalWntSignalingPathway <-"GO:0060070"
NonCanonicalWntSignalingPathway <- "GO:0035567"



