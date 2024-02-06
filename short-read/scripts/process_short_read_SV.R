

#########################################################
####### processing short read sequencing data ###########

rm(list = ls())

### load required libraries
library(dplyr)
library(ggplot2)
library(vcfR)
library(StructuralVariantAnnotation)
library(plyr)
library(UpSetR)



###########################
######### tsting ##########
###########################
### variants detected per SV caller
#all_calls<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/filtered-vcfs", pattern = "*.vcf", full.names = TRUE)


#vcfs_all<-lapply(all_calls,function(x){
#      aa<-VariantAnnotation::readVcf(x)
#      colo829_bpgr <- data.frame(c(breakpointRanges(aa),breakendRanges(aa)))
#
#})


#for (i in 1:length(vcfs_all)){
#    vcfs_all[[i]]<-cbind(vcfs_all[[i]],all_calls[i])
#    }

#vcfs_combined <- do.call("rbind", vcfs_all)


#names(vcfs_combined)[19]<-"path"
#ChrNames <- paste("chr",c(1:22,"X","Y") ,sep="")

#type_data_good<-vcfs_combined%>% 
#    mutate(sample=sub('.*/\\s*', '', gsub("/Users/shsoudi/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/filtered-vcfs/","",path)), )  %>% 
#    mutate(sample=gsub(".vcf","",sample)) %>% 
#    #filter((seqnames%in%ChrNames)) %>% 
#    distinct(sourceId , .keep_all= TRUE) %>% 
#    #filter(FILTER!="LowQual") %>% 
#    dplyr::select(!c(REF,ALT,path)) 
    

   

#vi_plot<-ggplot(type_data_good, aes(x = svtype, y = abs(svLen), fill = svtype)) +
#     geom_violin() +
#     scale_y_log10() +
#     geom_point(position = position_jitter(seed = 1, width = 0.2),alpha=0.01) +
#     theme(legend.position = "none") +
#    facet_wrap(~sample, nrow = 3)


#ggplot(delly, aes(x = svtype, y = abs(svLen))) +
#  geom_violin() +
#  scale_y_log10() +
#  geom_point(position = position_jitter(seed = 1, width = 0.2),alpha=0.1,aes(colour = svtype)) +
#  theme(legend.position = "none")   
    
##############################################################
####### Process individual vcf-bedpe files - Illumina ########
#############################################################

### load bedpe files (variants detected by each SV caller) 
unmerged_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_3_match",pattern = "*.bed", full.names= TRUE) 


unmerged_data<-lapply(unmerged_file, function(x){

    bed<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(unmerged_data)){
    unmerged_data[[i]]<-cbind(unmerged_data[[i]],unmerged_file[i])
    }

beds_illu<-do.call("rbind",unmerged_data)
names(beds_illu)[length(names(beds_illu))]<-"path"

beds_illu_good<-beds_illu%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100_3match.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    mutate(unique_id =paste(sample,V11, sep = "_")) %>%  ### V11 is the SV type
    mutate(SV_length =V2-V5) 
    
names(beds_illu_good)[names(beds_illu_good) == "V11"] <- "svtype"    

## plot vilin-dor plot  
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_3_match/plots/violin_SV_frequency_Illumina_3caller_intersection.pdf", width = 14, height = 18)
vi_Illuplot1<-ggplot(beds_illu_good, aes(x = svtype, y = abs(SV_length), fill = svtype)) +
     geom_violin() +
     scale_y_log10() +
     geom_point(position = position_jitter(seed = 1, width = 0.2),alpha=0.08) +
     theme(legend.position = "none") +
     xlab("SV type") + 
     ylab("Frequency") +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"),
     legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
     plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))) +
    facet_wrap(~sample, nrow = 3)
print(vi_Illuplot1)
dev.off()

##########################################################
####### Process survivor merged output - Illumina ########
##########################################################

### load "bedpe" files from merged survivor vcfs #####
bed_files<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_3_match",pattern = "*.bed", full.names= TRUE) 


beds<-lapply(bed_files, function(x){

    bed.data<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(beds)){
    beds[[i]]<-cbind(beds[[i]],bed_files[i])
    }

beds_all<-do.call("rbind",beds)
names(beds_all)[12]<-"path"

bedpe_good<-beds_all%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100_3match.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    mutate(unique_id =paste(sample,V11, sep = "_"))  ### V11 is the SV type
    


counts1 <- ddply(bedpe_good, .(bedpe_good$unique_id,bedpe_good$V11), nrow)    #### count number of SV per cell line
colnames(counts1)<-c("sampleid_SV","SV_type","SV_count") 
counts1$cell_line<-gsub("_.*$","",counts1$sampleid_SV)
     
### basic stack bar chart (SV type per frequency)
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_3_match/plots/stacked_barplot1_suvivor_3intersection_SV_frequency_Illumina.pdf", width = 10, height = 9)
plot_freq1<-ggplot(counts1,aes(x = SV_type, y = SV_count, fill = cell_line)) +  
#ggplot(counts1,aes(x = cell_line, y = SV_count, fill = SV_type)) 
     geom_bar(stat = "identity",position = "stack", width = 0.6) +
     coord_flip() +
     #scale_fill_brewer(palette = 12)
     scale_fill_manual(values = c("forestgreen", "brown2", "#3C518F")) +
     xlab("SV type") + 
     ylab("Frequency") +
     ggtitle("Frequency of SVs per cell line-Illumina") +
     geom_text(aes(label = SV_count),color="black",size=2,position = position_stack(vjust = .5)) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"),
     legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
     plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
print(plot_freq1)
dev.off()     



### basic stack bar chart (cell line per frequency)
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_3_match/plots/stacked_barplot2_suvivor_3intersection_frequency_Illumina.pdf", width = 10, height = 9)
plot_freq2<-ggplot(counts1,aes(x = cell_line, y = SV_count, fill = SV_type)) +  
#ggplot(counts1,aes(x = cell_line, y = SV_count, fill = SV_type)) 
     geom_bar(stat = "identity",position = "stack", width = 0.7) +
     coord_flip() +
     #scale_fill_brewer(palette = 12)
     scale_fill_manual(values = c("#CA463F", "#E2CA58", "#B0D094", "#67BBD2", "#3C518F")) +
     #scale_fill_brewer(palette = 12) +
     xlab("Sarcome cell line") + 
     ylab("Frequency") +
     ggtitle("Frequency of SVs per cell line-Illumina") +
     geom_text(aes(label = SV_count),color="black",size=3,position = position_stack(vjust = .5)) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"),
     legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
     plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
print(plot_freq2)
dev.off()      




    

#cars_by_cylinders_gears <- mtcars %>%
#  group_by(cyl, gear) %>%
#  summarise(count = n())

###################################
### survivir merged vcf file #####
##################################

#survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor", pattern = "*.vcf", full.names = TRUE)

#vcfs_survivor<-lapply(survivor,function(x){
#     #readVcf(x)
#     #info(readVcf(x))
##     data.frame(info(readVcf(x)))
#})


#for (i in 1:length(vcfs_survivor)){
#    vcfs_survivor[[i]]<-cbind(vcfs_survivor[[i]],survivor[i])
#   }

#vcfs_survivor_combined <- do.call("rbind", vcfs_survivor)

#names(vcfs_survivor_combined)[length(names(vcfs_survivor_combined))]<-"path"

#survivor_combined_good<-vcfs_survivor_combined%>% 
#    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100.vcf","",path)), )  %>% 
#    mutate(row_id=rownames(vcfs_survivor_combined)) %>% 
#    dplyr::select(!c(path)) 
    
#col_ids<-c("sample","SVTYPE","CHR2","row_id")    
#survivor_combined_good$uniq_id<- apply( survivor_combined_good[ , col_ids] , 1 , paste , collapse = "-" )


###### load bedpe files ####
beds_2<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor",pattern = "*.bed", full.names= TRUE) 


beds_surivor<-lapply(beds_2, function(x){

    bed.data<-read.delim(x, header = FALSE, sep = "\t")
})


for (i in 1:length(beds_surivor)){
    beds_surivor[[i]]<-cbind(beds_surivor[[i]],beds_2[i])
    }



beds_survivor_all<-do.call("rbind",beds_surivor)
names(beds_survivor_all)[length(names(beds_survivor_all))]<-"path"

bedpe_survivor_good<-beds_survivor_all%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    mutate(svLen=V2-V5)
    
col_ids<-c("sample","V11","V1","V7")    
bedpe_survivor_good$uniq_id<- apply( bedpe_survivor_good[ , col_ids] , 1 , paste , collapse = "-" )


### merge vcf and bed file


### plot violin 
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor/violin_suvivor_2_SV_frequency_Illumina.pdf", width = 10, height = 14)
plot_vi1<-vi_plot<-ggplot(bedpe_survivor_good, aes(x = V11, y = abs(svLen), fill = V11)) +
     geom_violin() +
     scale_y_log10() +
     geom_point(position = position_jitter(seed = 1, width = 0.2),alpha=0.06) +
     xlab("SV type") + 
     ylab("Frequency") +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"),
     legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
     plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10))) + 
    facet_wrap(~sample, nrow = 3)
print(plot_vi1)   
dev.off() 



#########################################
### upset plot "one match" survivor #####
#########################################
### merged vcf file 
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_1_match", pattern = "*.vcf", full.names = TRUE)

### lead info field
vcfs_survivor<-lapply(files_survivor,function(x){
     #readVcf(x)
     #info(readVcf(x))
     data.frame(info(readVcf(x)))
})


for (i in 1:length(vcfs_survivor)){
    vcfs_survivor[[i]]<-cbind(vcfs_survivor[[i]],files_survivor[i])
    }
type_data <- do.call("rbind", vcfs_survivor) 
rownames(type_data)<-NULL
names(type_data)[length(names(type_data))]<-"path"


type_data_good<-type_data%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100_1match.vcf","",path)), )  %>% 
    dplyr::select(-(path)) 



### load matrix files 
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_1_match", pattern = "*.txt", full.names = TRUE,recursive=TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
mat <- do.call("rbind", matrix_tables) 
names(mat)[length(names(mat))]<-"raw_id"


mat_good<-mat%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100_1match_matrix.txt","",raw_id)), )  %>% 
    dplyr::select(-(raw_id)) 


names(mat_good)[1:3]<-c("delly","manta","smoove")




vcf_matrix<-cbind(type_data_good,mat_good)
# columns to paste together
all_samples<-unique(vcf_matrix$sample)



for (jj in 1:length(all_samples)){


      vcf_matrix_focal<-vcf_matrix[vcf_matrix$sample ==all_samples[jj],]

      df_upset<-vcf_matrix_focal[,c("delly","manta","smoove")]

      pdf(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/upsetplot_SVcount_short_",all_samples[jj],".pdf",sep = ""),height= 8, width = 8)
      upset(df_upset, main.bar.color = "brown3",
      matrix.color="#299807", point.size=5,
      number.angles = 30,
      line.size = 1.5,
      text.scale = c(1.3, 1.3, 1),
      sets.bar.color=c("maroon","blue","orange"))
      dev.off()


} ### jj loop


############
#### END ###
############


###### related to long read
#### ** NOTE ** -> should move to long read-sequencing script (temporary here)

#matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/survivor/pacbio_survivor_1matched_SVLEN100_DIS100", pattern = "*.txt", full.names = TRUE)
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio2_patient", pattern = "*.txt", recursive = TRUE, full.names = TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
mat <- do.call("rbind", matrix_tables) 
names(mat)[length(names(mat))]<-"raw_id"


mat_good<-mat%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_merged_SVLEN100_DIS100_matrix.txt","",raw_id)), )  %>% 
    dplyr::select(-(raw_id)) 

#mat_good$sample <- gsub("1", "GCT",
#        gsub("2", "RD",
#        gsub("3", "SW", mat_good$sample)))
   

names(mat_good)[1:4]<-c("cuteSV","nanoSV","sniffles","svim")  
all_samples<-unique(mat_good$sample)


for (jj in 1:length(all_samples)){


      vcf_matrix_focal<-mat_good[mat_good$sample ==all_samples[jj],]

      df_upset<-vcf_matrix_focal[,c("cuteSV","nanoSV","sniffles","svim")]

      pdf(file = paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio2_patient/plots/upsetplot_SVcount_longread_",all_samples[jj],".pdf",sep = ""),height= 8, width = 8)
      upset(df_upset, main.bar.color = "brown3",
      matrix.color="#299807", point.size=5,
      number.angles = 10,
      line.size = 1.5,
      text.scale = c(1.3, 1.3, 1),
      sets.bar.color=c("maroon","blue","orange","#2c0107"))
      dev.off()


} ### jj loop




### load bedfiles 
### load "bedpe" files from merged survivor vcfs #####
bed_files<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio2_patient", pattern = "*.bed", recursive = TRUE, full.names = TRUE)


beds<-lapply(bed_files, function(x){

    bed.data<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(beds)){
    beds[[i]]<-cbind(beds[[i]],bed_files[i])
    }

beds_all<-do.call("rbind",beds)
names(beds_all)[length(names(beds_all))]<-"path"

bedpe_good<-beds_all%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_merged_SVLEN100_DIS100.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    #mutate(sample=gsub("1", "GCT",
        #gsub("2", "RD",
        #gsub("3", "SW", sample)))) %>% 
        mutate(unique_id =paste(sample,V11, sep = "_")) ### V11 is the SV type
    


counts1 <- ddply(bedpe_good, .(bedpe_good$unique_id,bedpe_good$V11), nrow)    #### count number of SV per cell line
colnames(counts1)<-c("sampleid_SV","SV_type","SV_count") 
counts1$cell_line<-gsub("_.*$","",counts1$sampleid_SV)
     
### basic stack bar chart (SV type per frequency)
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio2_patient/plots/stacked_barplot1_suvivor_3_intersection_SV_frequency_pacbio.pdf", width = 10, height = 9)
plot_freq1<-ggplot(counts1,aes(x = SV_type, y = SV_count, fill = cell_line)) +  
#ggplot(counts1,aes(x = cell_line, y = SV_count, fill = SV_type)) 
     geom_bar(stat = "identity",position = "stack", width = 0.6) +
     coord_flip() +
     #scale_fill_brewer(palette = 12)
     scale_fill_manual(values = c("forestgreen", "brown2", "#3C518F")) +
     xlab("SV type") + 
     ylab("Frequency") +
     ggtitle("Frequency of SVs per cell line-longread pacbio") +
     #geom_text(aes(label = SV_count),color="black",size=2,position = position_stack(vjust = .5)) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"),
     legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
     plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
print(plot_freq1)
dev.off()     



### basic stack bar chart (cell line per frequency)
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/long_reads/pacbio2_patient/plots/stacked_barplot2_suvivor_3caller_intersection_SV_frequency_pacbio.pdf", width = 10, height = 9)
plot_freq2<-ggplot(counts1,aes(x = cell_line, y = SV_count, fill = SV_type)) +  
#ggplot(counts1,aes(x = cell_line, y = SV_count, fill = SV_type)) 
     geom_bar(stat = "identity",position = "stack", width = 0.7) +
     coord_flip() +
     #scale_fill_brewer(palette = 12)
     scale_fill_manual(values = c("#CA463F", "#E2CA58", "#B0D094", "#67BBD2", "#3C518F")) +
     #scale_fill_brewer(palette = 12) +
     xlab("Sarcome cell line") + 
     ylab("Frequency") +
     ggtitle("Frequency of SVs per cell line-longread pacbio") +
     #geom_text(aes(label = SV_count),color="black",size=3,position = position_stack(vjust = .5)) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
     panel.background = element_blank(), axis.line = element_line(colour = "black"),
     legend.key = element_rect(fill = NA, colour = NA, linewidth = 0.25),axis.text=element_text(size=16),axis.title=element_text(size=18),
     plot.margin = margin(1, 1., 1, 1.5, "cm"),axis.title.y = element_text(margin = margin(t = 70, r = 20, b = 50, l = 10)))
print(plot_freq2)
dev.off()      




    