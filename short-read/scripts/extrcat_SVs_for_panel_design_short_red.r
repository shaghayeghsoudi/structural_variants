#rm(list = ls())

####################################################################################
### extrcat variants for the panel design (from vcf files genearted by survivor) ###
### works for short-read sequencing ###
####################################################################################

#library(webr)
library(dplyr)
#library(vcfR)
#library(tidyverse)
#library(StructuralVariantAnnotation)
library(VariantAnnotation)
#library(data.table)
#library(stringi)


### load vcf files
files_survivor<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines/panel/survivor_1matched_SVLEN100_DIS100/vcfs", pattern = "*_illumina_SVLEN100_DIS100_1caller.vcf", full.names = TRUE)

### lead rowranges field
#var_vcf<-function(x){
#     vcf<-readVcf(x)
#     #info(readVcf(x))
#     vcf_row<-data.frame(rowRanges(vcf))
#}
#
#lapply(files_survivor,var_vcf)

##############################
### extract "info" columns ###
##############################
vcfs_survivor<-lapply(files_survivor,function(x){
     vcf<-readVcf(x)
     #info(readVcf(x))
     vcf_row<-data.frame(rowRanges(vcf))
     #data.frame(info(readVcf(x)))
})

for (i in 1:length(vcfs_survivor)){
    vcfs_survivor[[i]]<-cbind(vcfs_survivor[[i]],files_survivor[i])
    }
type_data <- do.call("rbind", vcfs_survivor) 
names(type_data)[11]<-"path"

type_data_good<-type_data%>% 
    mutate(sample_info=sub('.*/\\s*', '', gsub("_illumina_SVLEN100_DIS100_1caller.vcf","",path)), ) %>%  
    mutate(sample_info=gsub(".vcf","",sample_info)) 

#######################
#### load vcf info #####
vcfs_info<-lapply(files_survivor,function(x){
     vcf<-readVcf(x)
     #info(readVcf(x))
     #vcf_row<-data.frame(rowRanges(vcf))
     data.frame(info(readVcf(x)))
})

#for (i in 1:length(vcfs_info)){
#    vcfs_info[[i]]<-cbind(vcfs_info[[i]],files_survivor[i])
#    }


type_data_info <- do.call("rbind", vcfs_info) 
rownames(type_data_info )<-NULL

type_data_inforow<-cbind(type_data_good,type_data_info)  ### final merged set
type_data_inforow<-type_data_inforow[,colnames(type_data_inforow)!="path"] %>% 
    mutate(sample_info=gsub("2callers","2caller",sample_info))



#########################
### load matrix files ###
#########################
matrix<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines/panel/survivor_1matched_SVLEN100_DIS100/matrix", pattern = "*.txt", full.names = TRUE)


matrix_tables<-lapply(matrix,function(x){
    read.csv(x, header = FALSE, strip.white=T, sep='')
})

for (i in 1:length(matrix_tables)){
    matrix_tables[[i]]<-cbind(matrix_tables[[i]],matrix[i])
    }
mat <- do.call("rbind", matrix_tables) 
names(mat)[5]<-"raw_id"

good_mat<-mat %>% mutate(sample_info=sub('.*/\\s*', '',gsub("_illumina_SVLEN100_DIS100_1caller_matrix.txt","",raw_id)))
#good_mat$sample_info<-gsub("overlapped_","",gsub("_matrix.txt","",good_mat$sample_info))

good_mat$total_intersection<-rowSums(good_mat[ , c(1:4)], na.rm=TRUE)

good_mat_fin<-good_mat%>%
        rename(V1="Delly" , V2="Manta" , V3= "Gridss",V4="Smoove") %>% 
        mutate(caller_count=rowSums(.[1:4])) %>% 
        mutate(sample_info=gsub("ways","caller",sample_info)) %>% 
        dplyr::select(-raw_id)
       


both<-mutate(type_data_inforow,good_mat_fin) 

chroms<-c(paste("chr",1:22, sep = ""),"chrX")
both<-both[(both$seqnames%in%chroms) & (both$CHR2%in%chroms),]
#both$caller_count_info<-sub('.*\\_', '', both$sample_info)

both$f<-substr(both$STRANDS, 1,1)
both$r<-substr(both$STRANDS, 2,2)


both$brachet_info<-ifelse(both$SVTYPE=="TRA",both$ALT,"NA")
both$ins_seq<-ifelse(both$SVTYPE=="INS",both$ALT,"NA")
selected_both<-both[,c("SVTYPE","seqnames","start","f","CHR2","END","r","Delly","Manta","Gridss","Smoove","sample_info")]


### extrcat bracet info and 
brac_final<-data.frame("brachet_info"=unlist(both$brachet_info))
brac_ins<-data.frame("ins_seq"=unlist(both$ins_seq))

both_final<-cbind(selected_both,brac_final,brac_ins)
write.table(both_final,file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines/panel/survivor_1matched_SVLEN100_DIS100,panel_shortread.txt",sep = "\t",quote = FALSE, col.names = TRUE, row.names = FALSE)

 ##############
 ##### END ####   


