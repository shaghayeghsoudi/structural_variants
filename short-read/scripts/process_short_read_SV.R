

#########################################################
####### processing short read sequencing data ###########

rm(list = ls())

### load required libraries
library(dplyr)
library(ggplot2)
library(vcfR)
library(StructuralVariantAnnotation)
library(plyr)



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
unmerged_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/filtered-bedfiles",pattern = "*.bed", full.names= TRUE) 


unmerged_data<-lapply(unmerged_file, function(x){

    bed<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(unmerged_data)){
    unmerged_data[[i]]<-cbind(unmerged_data[[i]],unmerged_file[i])
    }

beds_illu<-do.call("rbind",unmerged_data)
names(beds_illu)[length(names(beds_illu))]<-"path"

beds_illu_good<-beds_illu%>% 
    mutate(sample=sub('.*/\\s*', '', gsub(".filtered.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    mutate(unique_id =paste(sample,V11, sep = "_")) %>%  ### V11 is the SV type
    mutate(SV_length =V2-V5) 
    
names(beds_illu_good)[names(beds_illu_good) == "V11"] <- "svtype"    

## plot vilin-dor plot  
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor/violin_SV_frequency_Illumina.pdf", width = 14, height = 18)
vi_Illuplot1<-ggplot(beds_illu_good, aes(x = svtype, y = abs(SV_length), fill = svtype)) +
     geom_violin() +
     scale_y_log10() +
     geom_point(position = position_jitter(seed = 1, width = 0.2),alpha=0.01) +
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
bed_files<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor",pattern = "*.bed", full.names= TRUE) 


beds<-lapply(bed_files, function(x){

    bed.data<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(beds)){
    beds[[i]]<-cbind(beds[[i]],bed_files[i])
    }

beds_all<-do.call("rbind",beds)
names(beds_all)[12]<-"path"

bedpe_good<-beds_all%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_survivor_filt_merged_SVLEN100_DIS100.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    mutate(unique_id =paste(sample,V11, sep = "_"))  ### V11 is the SV type
    


counts1 <- ddply(bedpe_good, .(bedpe_good$unique_id,bedpe_good$V11), nrow)    #### count number of SV per cell line
colnames(counts1)<-c("sampleid_SV","SV_type","SV_count") 
counts1$cell_line<-gsub("_.*$","",counts1$sampleid_SV)
     
### basic stack bar chart (SV type per frequency)
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor/plots/stacked_barplot1_suvivor_2overlapped_SV_frequency_Illumina.pdf", width = 10, height = 9)
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
pdf("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor/plots/stacked_barplot2_suvivor_2_SV_frequency_Illumina.pdf", width = 10, height = 9)
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



#######################################
### upset plot one match survivor #####
#######################################


### load matrix files 
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_structural_variants/short-read/downstream2/vcfs/merged_survivor_1_match", pattern = "*.txt", full.names = TRUE,recursive=TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
good_mat <- do.call("rbind", matrix_tables) 
names(good_mat)[length(names(good_mat))]<-"raw_id"


good_mat$sample_id<-sub('.*/\\s*', '', good_mat$raw_id)
good_mat$sample_id<-gsub("survivor_merged_filtered_","",gsub("_overlapped_3ways_matrix.txt","",good_mat$sample_id))

good_mat_fin<-good_mat%>%
        separate(sample_id,c("sample","technology","minimum_SVLEN","distance"))


good_mat_fin$sample <- gsub("1", "GCT",
        gsub("2", "RD",
        gsub("3", "SW982", type_data_good$sample)))


names(good_mat_fin)[1:4]<-c("cuteSV","nanoSV","sniffles","svim")
good_mat_fin<-good_mat_fin[,-5]



vcf_matrix<-cbind(type_data_good,good_mat_fin)
# columns to paste together
all_samples<-unique(vcf_matrix$uniq_ID)


