

### find intersection between selected panel from longread and short read SV
rm(list = ls())

### load required libraries
library(dplyr)
library(ggplot2)
library(vcfR)
library(StructuralVariantAnnotation)
library(plyr)
library(UpSetR)
library(tidyr)

### load bedpe files (variants detected by each sjort read SV caller on resequenced cell lines) 
vcf_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read",pattern = "*.bedpe", full.names= TRUE) 

vcfs<-lapply(vcf_file, function(x){

    bed<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(vcfs)){
    vcfs[[i]]<-cbind(vcfs[[i]],vcf_file[i])
    }

beds_illu<-do.call("rbind",vcfs)
names(beds_illu)[length(names(beds_illu))]<-"path"

beds_illu_good<-beds_illu%>% 
    mutate(sample_caller=sub('.*/\\s*', '', gsub("_survivor.bedpe","",path)), )  %>% 
    dplyr::select(-c(path,V8)) %>% 
    mutate(sample=gsub("_.*$","",sample_caller)) %>% 
    dplyr::mutate(across(V2:V3,~.-1)) %>% 
    dplyr::mutate(across(V5:V6,~.-1)) 
    #mutate(unique_id =paste(sample,V11, sep = "_"))  ### V11 is the SV type
   
    
colnames(beds_illu_good) <- c("chrom1","start1","end1","chrom2","start2","end2","typeID","strand1","strand2","SVtype","sample_caller","sample") 
chroms_good<-c(paste("chr",seq(1:22),sep=""),"chrX","chrY")

beds_illu_good_ch<-beds_illu_good[(beds_illu_good$chrom1 %in% chroms_good) | (beds_illu_good$chrom2 %in% chroms_good),]


#########################
### load panel files ####
#########################
panel_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_panel",pattern = "*.bedpe", full.names= TRUE) 

panels<-lapply(panel_file, function(x){

    pa<-read.delim(x, header = TRUE, sep = "\t")
})

for (i in 1:length(panels)){
    panels[[i]]<-cbind(panels[[i]],panel_file[i])
    }

pan<-do.call("rbind",panels)
names(pan)[length(names(pan))]<-"path"

pan_good<-pan%>% 
    mutate(sample=sub('.*/\\s*', '', gsub("_only_finalprobebedfile.bedpe","",path)), )  %>% 
    dplyr::select(-(path)) %>% 
    dplyr::mutate(across(ChrPosA,~.+1)) %>% 
    dplyr::mutate(across(ChrPosB,~.+1)) %>% 
    mutate(SVtype = toupper(SVtype)) %>% 
    mutate(uniq_identifier=paste(sample,SVtype,ChrA,ChrPosA,ChrB,ChrPosB, sep = "_"))
    
    
pan_uniq<-pan_good %>% distinct(uniq_identifier, .keep_all= TRUE) 
    #mutate(by_sort=gsub("chr","",ChrA)) 


pan_good<-pan_good[order(as.numeric(pan_good$by_sort),decreasing = FALSE),]

pan_good_ch<-pan_good[(pan_good$ChrA%in% chroms_good | pan_good$ChrB%in% chroms_good),]   
    
samples<-unique(pan_good$sample)


for(ss in 1:length(samples)){

    focal_caller<-beds_illu_good_ch %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(chrom1,start1,start2,SVtype,sample_caller)
    #setDT(focal_caller_format)
    #setkey(focal_caller_format, chrom1,start1,start2)
    
    
    focal_panel<-pan_good_ch %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(ChrA,ChrPosA,ChrPosB,SVtype) 

    colnames(focal_panel)<-c("chrom1","start1","start2","SVtype")
    
    #rename_with(pan_good_ch, recode, ChrA = "chrom1",ChrPosA ="start1" ,  ChrPosB="start2" , SVtype= "SVtype")  ### fix renames


    chroms<-unique(focal_panel$chrom1)
    out_res<-NULL
    for(ii in 1:length(chroms)){

        focal_caller_chrom<-focal_caller[focal_caller$chrom1 ==chroms[ii], ] 
        focal_caller_chrom<-focal_caller[focal_caller$SVtype != "TRA",]       
        
        
        focal_panel_chrom<-focal_panel[focal_panel$chrom1 ==chroms[ii],]
        focal_panel_chrom<-focal_panel[focal_panel$SVtype != "TRA",]    


        setDT(focal_caller_chrom)
        setDT(focal_panel_chrom)

        
        
        setkey(focal_caller_chrom, chrom1,start1,start2)
        setkey(focal_panel_chrom, chrom1,start1,start2)

        overlaps<-data.frame(foverlaps(focal_caller_chrom, focal_panel_chrom, type="any"))
        overlaps_matchSV<-overlaps[overlaps$SVtype == overlaps$i.SVtype,] %>% 
           drop_na()
       out_res<-rbind(overlaps_matchSV,out_res) 




    }





    
    setDT(focal_panel_format)
    setkey(focal_panel_format, chrom1,start1,start2)


    overlaps <- foverlaps(focal_caller_format, focal_panel_format_sorted, type="any")





}

### from here
### create matrix
results_out <- array (NA, c((nrow (input_good) * (end1 - start1 + 1)),7))
count <- 0

system.time(
#loop through focal snp
for (i in start1:end1){
	
#	loop through all SNPs
	for (j in 1:nrow (input_good)){
		count <- count + 1
		 results_out[count,7] <- cor (input_good[i,], input_good[j,], use = "pairwise.complete.obs")
		results_out [count,1] <- scafpos_good[i,1]
		results_out [count,2] <- scafpos_good[i,2]
		#results_out [count,3] <- scafpos_good[i,3]
		results_out [count,4] <- scafpos_good[j,1]
		results_out [count,5] <- scafpos_good[j,2]
		#results_out [count,6] <- scafpos_good[j,3]
	
	}
}
)


outname <- paste ("output", start1,"_",end1,".txt",sep = "")
write.table (results_out, outname, col.names = F, row.names = F, quote = F)
