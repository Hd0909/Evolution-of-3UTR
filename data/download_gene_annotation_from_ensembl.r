#setwd("/media/huangdan/hardisk0/HD/my_data/HD_evolution/final_program_UTR9_21/")
#### It may take some time to download the data from ensembl database
library(biomaRt)
library(DescTools)
source('functions/functions_evolution.r')
#file contains the TSGs , Oncogenes , Non-cancer genes  list we used in our previous paper[2].The TSG and oncogene lists were obtained from the TSGene2.0[3] database and Vogelstein's review[4], respectively, as we previously described[3,5-6]. Genes not in these lists were considered as non-cancer genes. 
#[2],Wang, X. et al. Oncogenes expand during evolution to withstand somatic amplification. Ann. Oncol. 29, 2254–2260 (2018).
#[3],Zhao, M., Kim, P., Mitra, R., Zhao, J. & Zhao, Z. TSGene 2.0: An updated literature-based knowledgebase for Tumor Suppressor Genes. Nucleic Acids Res. 44, D1023–D1031 (2016).
#[4],Vogelstein, B. et al. Cancer genome landscapes. Science (80-. ). 340, 1546–1558 (2013).
#[5],Wu, W. K. K. et al. Oncogenes without a neighboring tumor-suppressor gene are more prone to amplification. Mol. Biol. Evol. 34, 903–907 (2017).
#[6],Huang, D. et al. Multi-omic analysis suggests tumor suppressor genes evolved specific promoter features to optimize cancer resistance. Brief. Bioinform. (2021) doi:10.1093/bib/bbab040.
#
all_genes <- read.csv("./data/all_genes.csv", stringsAsFactors=FALSE)
rownames(all_genes)<-all_genes[,1]
print(dim(all_genes))
print(head(all_genes))

ensembl = useEnsembl(biomart="ensembl","hsapiens_gene_ensembl",
                     host="asia.ensembl.org", version = 97, verbose = TRUE)
attributes.list<-listAttributes(ensembl)
att=c('ensembl_gene_id',"ensembl_transcript_id","ensembl_peptide_id","transcript_biotype","transcript_start","transcript_end","transcript_length","transcript_appris",'hgnc_symbol','chromosome_name','ensembl_exon_id', 'start_position', 'end_position', 'strand')

att2 <- attributes.list[ grep("aralog",attributes.list$description), 1]
# get the gene information in our gene list
get_gene_pos=getBM(attributes = att, mart = ensembl, 
                   filters=c('ensembl_gene_id'),
                   values=list(all_genes[,1]), uniqueRows = TRUE)

get_gene_pos_paralog=getBM(attributes = c(att[1],att2), mart = ensembl, 
                           filters=c('ensembl_gene_id'),
                           values=list(all_genes[,1]), uniqueRows = TRUE)
att3=c('ensembl_gene_id',"ensembl_transcript_id","ensembl_peptide_id",'5_utr_start','5_utr_end','3_utr_start','3_utr_end')
get_gene_pos_utr=getBM(attributes = att3, mart = ensembl, 
                   filters=c('ensembl_gene_id'),
                   values=list(all_genes[,1]), uniqueRows = TRUE)

get_gene_pos=get_gene_pos[which(get_gene_pos$chromosome_name %in% c(1:22,"X","Y")),]
get_gene_pos<-merge(get_gene_pos,get_gene_pos_utr)
print(head(get_gene_pos))
print(head(get_gene_pos_paralog))
## only selected protein-coding genes
get_gene_pos=get_gene_pos[which(get_gene_pos$transcript_biotype=="protein_coding"),]
get_gene_pos_paralog=get_gene_pos_paralog[which(get_gene_pos_paralog$ensembl_gene_id %in% get_gene_pos$ensembl_gene_id),]
all_data_genes<-list(get_gene_pos=get_gene_pos,gene_paralog=get_gene_pos_paralog)


########### file contains  age idntification of human genes and it was obtained from the supplementary table 2 of the Yin et al.'s paper[1]
########Yin, H., Wang, G., Ma, L., Yi, S. V. & Zhang, Z. What Signatures dominantly associate with gene age? Genome Biol. Evol. 8, 3083–3089 (2016).
gene_age <- read.delim("./data/gene_age.txt", stringsAsFactors=FALSE)
rownames(gene_age)<-gene_age$ID
gene_age_infor <- read.delim("./data/gene_age_infor.txt")
rownames(gene_age_infor)<-gene_age_infor$Age.class

print(str(gene_age))
print(summary(gene_age))

### data obtained from ensembl --> step1_download_from_ensembl_primary_data.r
print(str(all_data_genes$get_gene_pos))
print(summary(all_data_genes$get_gene_pos))

print(str(all_data_genes$gene_paralog))
print(summary(all_data_genes$gene_paralog))

delete_cols=which(colnames(all_data_genes$get_gene_pos) %in% c("ensembl_transcript_id","ensembl_peptide_id","transcript_biotype","transcript_start","transcript_end","ensembl_exon_id","transcript_appris","transcript_length",'5_utr_start','5_utr_end','3_utr_start','3_utr_end'))
genes_infor<-unique(all_data_genes$get_gene_pos[,-c(delete_cols)])
# checked this two genes are renamed or deleted 
genes_infor=genes_infor[(!genes_infor$hgnc_symbol%in% c("CCL3L3","PLEKHG7")) & genes_infor$hgnc_symbol!="",]
#3 TSG did not have hgnc symbol  AC016586.1 ENSG00000077463

####################################### May not use the transcript count in this project

genes_infor$all_trans_count<-0
genes_infor$pro_count<-0
genes_infor$non_pro_trans_count<-0


genes_infor$Age_class=gene_age[genes_infor[,1],]$Age_class
temp_age=rep("1_13",length(genes_infor$Age_class))
temp_age[which(genes_infor$Age_class>13 & genes_infor$Age_class<20)]="13_20"
temp_age[which(genes_infor$Age_class>=20 & genes_infor$Age_class<23)]="20_23"
temp_age[which(genes_infor$Age_class>=23 )]="23+"
temp_age[is.na(genes_infor$Age_class)]=NA
genes_infor$age_type=temp_age

gene_paralog=all_data_genes$gene_paralog
gene_paralog=gene_paralog[gene_paralog$hsapiens_paralog_ensembl_gene!="",]
gene_paralog_count=table(gene_paralog$ensembl_gene_id)

genes_infor$gene_paralog_count=as.numeric(gene_paralog_count[genes_infor$ensembl_gene_id])
genes_infor$Gene.Type=all_genes[genes_infor$ensembl_gene_id,]$Gene.Type

############################################################################################### get main transcript
get_gene_pos_temp<-unique(all_data_genes$get_gene_pos[,-(which(colnames(all_data_genes$get_gene_pos)=="ensembl_exon_id"))])
get_gene_pos_temp$Gene.Type=all_genes[get_gene_pos_temp$ensembl_gene_id,]$Gene.Type# ensembl ID of "CCL3L3","PLEKHG7" is out of date.Those genes have new names now
get_gene_pos_temp=get_gene_pos_temp[(!get_gene_pos_temp$hgnc_symbol%in% c("CCL3L3","PLEKHG7")) & get_gene_pos_temp$hgnc_symbol!="",]

get_gene_pos_temp$APPRIS.annotation<-get_gene_pos_temp$transcript_appris
get_gene_pos_temp$trans_length=get_gene_pos_temp$transcript_length
##### get the main transcript ; first select genes with APPRIS princinple 1 transcript, if no APPRIS annotation, select the transcript with max length. if two transcript have the same length , get the first one
genes_main_transcript<-get_main_transcript(get_gene_pos_temp)
rownames(genes_main_transcript)<-genes_main_transcript$ensembl_gene_id
genes_infor$genes_main_transcript<-genes_main_transcript[genes_infor$ensembl_gene_id,]$ensembl_transcript_id

genes_infor$Origin.time..million.years.ago.._original=gene_age_infor[genes_infor$ensembl_gene_id,"Origin.time..million.years.ago.."]

get_gene_pos=get_gene_pos[which(get_gene_pos$chromosome_name %in% c(1:22,"X","Y")),]

print(head(get_gene_pos))
print(head(get_gene_pos_paralog))

temp_gene_pos=get_gene_pos_utr
temp_gene_pos$len_5UTR=abs(get_gene_pos_utr$`5_utr_end`- get_gene_pos_utr$`5_utr_start`+1)
temp_gene_pos$len_3UTR=abs(get_gene_pos_utr$`3_utr_end`- get_gene_pos_utr$`3_utr_start`+1)

len_5UTR=aggregate(temp_gene_pos$len_5UTR,list(temp_gene_pos$ensembl_transcript_id),function(x){sum(x,na.rm=T)})
rownames(len_5UTR)<-len_5UTR$Group.1
len_3UTR=aggregate(temp_gene_pos$len_3UTR,list(temp_gene_pos$ensembl_transcript_id),function(x){sum(x,na.rm=T)})
rownames(len_3UTR)<-len_3UTR$Group.1

get_gene_pos$len_5UTR<-len_5UTR[get_gene_pos$ensembl_transcript_id,2]
get_gene_pos$len_3UTR<-len_3UTR[get_gene_pos$ensembl_transcript_id,2]
get_gene_pos$len_5UTR[get_gene_pos$len_5UTR==0]=NA
get_gene_pos$len_3UTR[get_gene_pos$len_3UTR==0]=NA
all_data_genes<-list(get_gene_pos=get_gene_pos,gene_paralog=get_gene_pos_paralog)

temp_x=unique(all_data_genes$get_gene_pos[,c("len_5UTR","len_3UTR","ensembl_transcript_id","ensembl_peptide_id")])
rownames(temp_x)<-temp_x$ensembl_transcript_id
genes_infor$hg38_5UTR_length<-temp_x[genes_infor$genes_main_transcript,]$len_5UTR
genes_infor$hg38_3UTR_length<-temp_x[genes_infor$genes_main_transcript,]$len_3UTR
genes_infor$ensembl_peptide_id<-temp_x[genes_infor$genes_main_transcript,]$ensembl_peptide_id
genes.length<-unique(all_data_genes$get_gene_pos[,c("ensembl_transcript_id","transcript_length")])
rownames(genes.length)<-genes.length[,1]
genes_infor$genes_trans_length<-genes.length[genes_infor$genes_main_transcript,2]

genes_infor$ratio_3UTR_5UTR=genes_infor$hg38_3UTR_length/genes_infor$hg38_5UTR_length
genes_infor$cds_length<-genes_infor$genes_trans_length-genes_infor$hg38_3UTR_length-genes_infor$hg38_5UTR_length

genes_infor$hg38_3UTR_length_log<-log10(genes_infor$hg38_3UTR_length)
genes_infor$hg38_5UTR_length_log<-log10(genes_infor$hg38_5UTR_length)
genes_infor$cds_length_log<-log10(genes_infor$cds_length)
genes_infor$hg38_3utr_ratio_all<-genes_infor$hg38_3UTR_length/genes_infor$genes_trans_length
genes_infor$genes_trans_length_log<-log10(genes_infor$genes_trans_length)
genes_infor$cds_ratio_all<-genes_infor$cds_length/genes_infor$genes_trans_length
genes_infor$hg38_3utr_ratio_non5UTR<-genes_infor$hg38_3UTR_length/(genes_infor$genes_trans_length-genes_infor$hg38_5UTR_length)
genes_infor$hg38_3utr_to_cds<-log10(genes_infor$hg38_3UTR_length/genes_infor$cds_length)
genes_infor$hg38_5utr_to_cds<-log10(genes_infor$hg38_5UTR_length/genes_infor$cds_length)
genes_infor$log_ratio_3UTR_5UTR<-log10(genes_infor$ratio_3UTR_5UTR)

genes_infor<-genes_infor %>% mutate(len_3UTR_log_bins=cut(log10(genes_infor$hg38_3UTR_length), breaks=10))
genes_infor$len_3UTR_log_bins=as.character(genes_infor$len_3UTR_log_bins)
genes_infor$len_3UTR_log_bins[which(genes_infor$len_3UTR_log_bins=="(-0.00452,0.452]")]="(0,0.452]"

genes_infor<-genes_infor %>% mutate(genes_trans_length_log=cut(genes_infor$genes_trans_length_log, breaks=10))
genes_infor<-genes_infor %>% mutate(ratio_all_3UTR_log_bins=cut(genes_infor$hg38_3utr_ratio_all, breaks=10))
genes_infor$ratio_all_3UTR_log_bins=as.character(genes_infor$ratio_all_3UTR_log_bins)
genes_infor$ratio_all_3UTR_log_bins[which(genes_infor$ratio_all_3UTR_log_bins=="(-0.000485,0.1]")]="(0,0.1] "



print("get the miRNA binding site in 3UTR from miRWalk database")
##########################################################################################
ensembl_NCBI_transcript_hg38 <- read.delim("./data/ensembl_NCBI_transcript_hg38.txt")
hsa_miRWalk_3UTR <- read.delim("./data/hsa_miRWalk_3UTR.txt")
hsa_miR_3UTR=unique(hsa_miRWalk_3UTR[,c(1:4)])
hsa_mir_3utr_1<-unique(hsa_miR_3UTR[,c(1,2)])
mirna_count_3utr<-as.data.frame(table(hsa_mir_3utr_1$mRNA))
rownames(mirna_count_3utr)<-mirna_count_3utr$Var1

hsa_mir_3utr_2<-unique(hsa_miR_3UTR[,c(2,4)])
mirna_binding_3utr<-as.data.frame(table(hsa_mir_3utr_2$mRNA))
rownames(mirna_binding_3utr)<-mirna_binding_3utr$Var1


ncbi_ens<-ensembl_NCBI_transcript_hg38[ensembl_NCBI_transcript_hg38$RefSeq.mRNA.ID !="" & ensembl_NCBI_transcript_hg38$Transcript.stable.ID %in% genes_infor$genes_main_transcript,]

ncbi_ens$mirna_count_3utr=mirna_count_3utr[ncbi_ens$RefSeq.mRNA.ID,2]
ncbi_ens$mirna_binding_3utr=mirna_binding_3utr[ncbi_ens$RefSeq.mRNA.ID,2]

mirna_count_ens=aggregate(ncbi_ens[,c("mirna_count_3utr","mirna_binding_3utr")],list(ncbi_ens$Gene.stable.ID),max)
rownames(mirna_count_ens)<-mirna_count_ens$Group.1
genes_infor$mirna_3UTR_count=mirna_count_ens[genes_infor$ensembl_gene_id,"mirna_count_3utr"]
genes_infor$mirna_binding_3utr<-mirna_count_ens[genes_infor$ensembl_gene_id,"mirna_binding_3utr"]



###################################################################################
print("Get RBP binding sites in 3UTR from beRBP database") 
#RBP-G_targets_26RBPs$ grep -H ">" *.fa > all_RBP_with_RBPnames.txt

all_utr_RBPs <- unique(read.csv("./data/beRBP-G_targets_26RBPs/all_RBP_with_RBPnames.txt", header=FALSE, sep=";"))
x1=strsplit(all_utr_RBPs[,1],"_",fixed=T)
x_1=do.call(rbind.data.frame,x1)
colnames(x_1)<-c("gene","motif")
x_1$position=all_utr_RBPs[,2]

temp=do.call(rbind,strsplit(x_1$gene,"\\.|>"))
x_1$gene=temp[,5]
x_1$RBP=temp[,2]


knownToEnsembl_19 <- read.delim("./data/knownToEnsembl.txt", header=FALSE)

x3=do.call(rbind.data.frame,strsplit(knownToEnsembl_19$V1,".",fixed=T))
knownToEnsembl_19$V1<-x3[,1]
rownames(knownToEnsembl_19)<-knownToEnsembl_19$V1
x_1$ens<-knownToEnsembl_19[x_1$gene,2]
gene_rbp=unique(x_1[!is.na(x_1$ens),c("position","ens","RBP")])
gene_rbp_support<-gene_rbp[grepl("supported_by",gene_rbp[,1]),]
rbp_site_binding2=as.data.frame(table(gene_rbp$ens))
rownames(rbp_site_binding2)<-rbp_site_binding2$Var1

gene_rbp_uni=unique(x_1[!is.na(x_1$ens),c("ens","RBP")])
rbp_type_count<-as.data.frame(table(gene_rbp_uni$ens))
rownames(rbp_type_count)<-rbp_type_count$Var1



genes_infor$rbp_site_binding2<-rbp_site_binding2[genes_infor$genes_main_transcript,2]


genes_infor$rbp_type_count<-rbp_type_count[genes_infor$genes_main_transcript,2]
genes_infor$rbp_site_binding2<-rbp_site_binding2[genes_infor$genes_main_transcript,2]
genes_infor$rbp_site_binding_density2<-log10(genes_infor$rbp_site_binding2/genes_infor$hg38_3UTR_length)

rbp_site_binding2_support=as.data.frame(table(gene_rbp_support$ens))
rownames(rbp_site_binding2_support)<-rbp_site_binding2_support$Var1

genes_infor$rbp_site_binding2_support<-rbp_site_binding2_support[genes_infor$genes_main_transcript,2]
genes_infor$rbp_site_binding_density2_support<-log10(genes_infor$rbp_site_binding2_support/genes_infor$hg38_3UTR_length)


genes_infor$mirna_3UTR_count_density=log10(genes_infor$mirna_3UTR_count/genes_infor$hg38_3UTR_length)
genes_infor$mirna_binding_3utr_density=log10(genes_infor$mirna_binding_3utr/genes_infor$hg38_3UTR_length)


genes_infor$hg38_5utr_ratio_all<-genes_infor$hg38_5UTR_length/genes_infor$genes_trans_length
genes_infor$hg38_3utr_ratio_all<-genes_infor$hg38_3UTR_length/genes_infor$genes_trans_length

genes_infor$hg38_3utr_ratio_cds<-genes_infor$hg38_3UTR_length/genes_infor$cds_length


print("genes_name")
rownames(genes_infor)<-genes_infor[,1]


save(genes_infor,all_data_genes,file=paste("./data/genes_infor_downloaded_from_ensembl97.rdata",sep=""))
##### print("other species")
####################################################### other species
path <- "./data/Species_we_selected_V2.csv"
Data<-read.csv(path)
Species.select<- Data$Species

Species.select<-Species.select[Species.select!=""]

ensembl_version=97
ensembl = useEnsembl(biomart="ensembl",
                     host="asia.ensembl.org", version = ensembl_version, verbose = TRUE)
en_Datasets<-listDatasets(ensembl)

ex <- paste(Species.select, collapse = " genes|^")
en_Datasets.selected<-en_Datasets[grep(ex, en_Datasets$description),]

print(ensembl_version)
for (i in 1:dim(en_Datasets.selected)[1])
{
  print(en_Datasets.selected[i,])
  ensembl_used = useEnsembl(biomart="ensembl",
                            dataset = as.character(en_Datasets.selected$dataset[i]),
                            host="asia.ensembl.org", version = ensembl_version, verbose = FALSE)
  attributes.list<-listAttributes(ensembl_used)
  print(en_Datasets[en_Datasets$dataset ==en_Datasets.selected[i,1],])
  a<-subset(attributes.list, page =="feature_page")
  att=c('ensembl_gene_id',"ensembl_transcript_id",'chromosome_name', "transcript_biotype",'start_position', 'end_position', 'strand','5_utr_start','5_utr_end','3_utr_start','3_utr_end')
  ## some species did not have transcript_appris
  att0=c("ensembl_transcript_id","transcript_start","transcript_end","transcript_appris","transcript_length")
  att0=att0[att0 %in% ensembl_used@attributes$name]
  fname = paste("UTR_Data_of_Feature_",
                strsplit( en_Datasets.selected$description[i], " genes ")[[1]][1], ".csv",sep = "")
  fpath = "./data/"
  fpath = paste(fpath,fname,sep = "")
  Data.getted = getBM(attributes = att, mart = ensembl_used, bmHeader = TRUE, uniqueRows = TRUE)
  
  Data.getted0 = getBM(attributes = att0, mart = ensembl_used, bmHeader = TRUE, uniqueRows = TRUE)
  re=merge(Data.getted,Data.getted0,by="Transcript stable ID")
  write.csv(re,fpath)
}




ensembl_used = useEnsembl(biomart="ensembl",
                          dataset = "hsapiens_gene_ensembl",
                          host="asia.ensembl.org", version = ensembl_version, verbose = TRUE)

attributes.list<-listAttributes(ensembl_used)
a<-subset(attributes.list, page =="homologs")
for (i in 1:length(Species.select))
{
  print(Species.select[i])
  descriptions <- c("gene stable ID", "gene name", "chromosome/scaffold name", "homology type")
  dscip_02 <- paste(Species.select[i], descriptions)
  att_02 <- a[a$description %in% dscip_02,1]
  att_01 <- c("ensembl_gene_id","external_gene_name",
              "chromosome_name")
  att <- c(att_01, att_02)
  fname = paste("Converter_Human_to_", att_02[1], "_V2.csv",sep = "")
  fpath = "./data/"
  fpath = paste(fpath,fname,sep = "")
  Data.getted = getBM(attributes = att, mart = ensembl_used, bmHeader = TRUE, uniqueRows = TRUE)
  write.csv(Data.getted,fpath)
}


dataset.name<-data.frame(name1=gsub("gene_ensembl","homolog_ensembl_gene",en_Datasets.selected$dataset), name2=do.call(rbind,strsplit( en_Datasets.selected$description, " genes "))[,1])


######################333 get the position of orhologs of other species
get_gene_pos_species_utr<-function(species_x,sp,genes_infor){
  Data_of_Feature <- read.csv(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/data/UTR_Data_of_Feature_",species_x,".csv",sep=""))
  ortholog <- read.csv(paste("/media/huangdan/hardisk0/HD/my_data/HD_evolution/data/Converter_Human_to_", sp, "_V2.csv",sep=""))
  ortholog=ortholog[!ortholog[,5]=="",]
  ortholog$Gene.Type=genes_infor[ortholog$Gene.stable.ID,]$Gene.Type
  TSG_non_cancer<-ortholog[which(ortholog$Gene.Type!="Oncogene"),]
  non_cancer_data=ortholog[which(ortholog$Gene.Type=="Non-Cancer"),]
  ## gene-noncancer pairs where gene-TSG exist at the same time
  non_cancer_dup=non_cancer_data[non_cancer_data[,5] %in% ortholog[ortholog$Gene.Type=="TSG",5],]
  
  #### delete the gene-noncancer pairs if gene-TSG exist
  unique_orthologs_select=unique(TSG_non_cancer[! paste(TSG_non_cancer[,2],TSG_non_cancer[,5],sep="_") %in% paste(non_cancer_dup[,2],non_cancer_dup[,5],sep="_"),c(5,9)])
  print("thus I keeep those genes and regard them as TSG")
  rownames(unique_orthologs_select)<-as.character(unique_orthologs_select[,1])
  if(!("APPRIS.annotation" %in% colnames(Data_of_Feature))){
    Data_of_Feature$APPRIS.annotation=rep("",length(Data_of_Feature$X))
  }
  get_gene_pos_temp<-Data_of_Feature[,c("Chromosome.scaffold.name","Transcript.stable.ID","Gene.stable.ID","Transcript.start..bp.","Transcript.end..bp.",                                "Strand","APPRIS.annotation","Transcript.length..including.UTRs.and.CDS.","Transcript.type")]
  colnames(get_gene_pos_temp)<-c("chromosome_name","ensembl_transcript_id","ensembl_gene_id","transcript_start","transcript_end","strand",
                                 "APPRIS.annotation","trans_length", "transcript_biotype")
  get_gene_pos_temp=get_gene_pos_temp[get_gene_pos_temp$transcript_biotype=="protein_coding",]
  genes_main_transcript<-get_main_transcript(get_gene_pos_temp)
  
  gene_pos=unique(genes_main_transcript[,c("ensembl_transcript_id","chromosome_name",                                   "transcript_start","transcript_end","strand","ensembl_gene_id","trans_length")])
  colnames(gene_pos)<-c("trans_ID","seqnames","start","end","strand","gene_id","trans_length")
  gene_pos$strand[gene_pos$strand=="1"]="+"
  gene_pos$strand[gene_pos$strand=="-1"]="-"
  gene_pos$seqnames[gene_pos$seqnames=="MT"]="M"
  gene_pos$seqnames=paste("chr",gene_pos[,2],sep="")
  # DELETE GENES WITH NO PROMOTER START =1 chrM  -999 100  1100      + ENSMUST00000082387 ENSMUSG00000064336
  if(species_x!="Human"){
    gene_pos<- gene_pos[gene_pos$gene_id %in% ortholog[,5],]}
  #gene_pos_data<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=T)
  Data_of_Feature$len_5UTR=abs(Data_of_Feature$X5..UTR.end- Data_of_Feature$X5..UTR.start+1)
  Data_of_Feature$len_3UTR=abs(Data_of_Feature$X3..UTR.end- Data_of_Feature$X3..UTR.start+1)
  len_5UTR=aggregate(Data_of_Feature$len_5UTR,list(Data_of_Feature$Transcript.stable.ID),function(x){sum(x,na.rm=T)})
  rownames(len_5UTR)<-len_5UTR$Group.1
  len_3UTR=aggregate(Data_of_Feature$len_3UTR,list(Data_of_Feature$Transcript.stable.ID),function(x){sum(x,na.rm=T)})
  rownames(len_3UTR)<-len_3UTR$Group.1
  
  print(head(gene_pos))
  result=data.frame(Group.1=gene_pos$gene_id,trans_length=gene_pos$trans_length, len_5UTR=len_5UTR[gene_pos$trans_ID,2],len_3UTR=len_3UTR[gene_pos$trans_ID,2],Gene.Type=unique_orthologs_select[gene_pos$gene_id,"Gene.Type"],species=species_x)
  return(result) 
}
result.utr.others<-c()
for(i in 1:length(dataset.name[,1])){
  print(dataset.name[i,1])
  if(dataset.name[i,2]!="Human"){
    temp<-get_gene_pos_species_utr(dataset.name[i,2],dataset.name[i,1],genes_infor)
    result.utr.others<-rbind(result.utr.others,temp)
  }
}

save(result.utr.others,file="./data/other_species_utr_more.rdata")
















