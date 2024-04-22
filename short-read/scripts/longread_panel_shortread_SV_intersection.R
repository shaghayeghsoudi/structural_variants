

### find intersection between selected panel from longread and short read SV
rm(list = ls())

### load required libraries
library(dplyr)
#library(ggplot2)
#library(vcfR)
#library(StructuralVariantAnnotation)
#library(plyr)
#library(UpSetR)
library(tidyr)
library(data.table)



############################################################
### load panel files obtained from longread resequencing####
############################################################

panel_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_panel",pattern = "*.bedpe", full.names= TRUE) 

## alternative
#my_fa<-function(x){  ## write the function first
#    read.delim(x, header = TRUE, sep = "\t")
#}
#panels<-lapply(panel_file, my_fa)  ## run in in lappy

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
    dplyr::mutate(across(ChrPosA,~.+1),across(ChrPosB,~.+1)) %>% 
    mutate(SVtype = toupper(SVtype),uniq_identifier=paste(sample,SVtype,ChrA,ChrPosA,ChrB,ChrPosB, sep = "_")) %>% 
    distinct(uniq_identifier, .keep_all= TRUE) 

#pan_uniq<-pan_uniq[order(as.numeric(pan_uniq$by_sort),decreasing = FALSE),] 

pan_chroms<-unique(c(pan_uniq$ChrA ,pan_uniq$ChrB))  ### unique chromosmes present in the pannel



###################################################################################################
### load bedpe files (variants detected by each sjort read SV caller on resequenced cell lines) ###
###################################################################################################

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
   
    
colnames(beds_illu_good) <- c("chromA","startA","endA","chromB","startB","endB","typeID","strand1","strand2","SVtype","sample_caller","sample") 
#chroms_good<-c(paste("chr",seq(1:22),sep=""),"chrX","chrY")

beds_illu_good_ch<-beds_illu_good[(beds_illu_good$chromA %in% pan_chroms) & (beds_illu_good$chromB %in% pan_chroms),]  ### remove unplaced scaffolds


### find overlaps     
samples<-unique(pan_uniq$sample)  ### what are the samples in the panel 


for(ss in 1:length(samples)){
  
    ## SV callers from resequenced 
    focal_caller<-beds_illu_good_ch %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(chromA,startA,startB,SVtype,sample_caller)
    #setDT(focal_caller_format)
    #setkey(focal_caller_format, chrom1,start1,start2)
    
    ## longread panel
    focal_panel<-pan_uniq %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(ChrA,ChrPosA,ChrPosB,SVtype) 

    colnames(focal_panel)<-c("chromA","startA","startB","SVtype")
    
    #rename_with(pan_good_ch, recode, ChrA = "chrom1",ChrPosA ="start1" ,  ChrPosB="start2" , SVtype= "SVtype")  ### fix renames


    chroms<-unique(focal_panel$chromA)
    
    out_res_chrom<-NULL
    for(ii in 1:length(chroms)){
        #for(ii in 1:20){

  
        ## SV caller
        focal_caller_chrom<-focal_caller %>%
         filter(chromA ==chroms[ii] & SVtype != "TRA") %>%  ### tTO FIX: emporarily skip TRA
         filter(!(startB < startA))              ### tTO FIX: emporarily skip TRA
         
           
        
        ## panel
        focal_panel_chrom<-focal_panel %>% 
         filter(chromA ==chroms[ii] & SVtype != "TRA")%>%  ### tTO FIX: emporarily skip TRA
         filter(!(startB < startA))

        setDT(focal_caller_chrom)
        setDT(focal_panel_chrom)

        
        
        setkey(focal_caller_chrom, chromA,startA,startB)
        setkey(focal_panel_chrom, chromA,startA,startB)

        overlaps<-data.frame(foverlaps(focal_caller_chrom, focal_panel_chrom, type="any"))
        overlaps_matchSV<-overlaps[overlaps$SVtype == overlaps$i.SVtype,] %>% 
           drop_na()

        perfect_match<-overlaps_matchSV[overlaps_matchSV$startA==overlaps_matchSV$i.startA & overlaps_matchSV$startB==overlaps_matchSV$i.startB,]   


           ### rbind chroms
       if (nrow(perfect_match)>0) {

            out_res_chrom<-rbind(perfect_match,out_res_chrom) 


       } ## if loop 


    }  ### chrom loop 

}    


write.table(out_res_chrom, file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/outputs/out_res_longread_pannel_resequenced_short_read_SW_other3.table", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)

#########################
###### just gridds #######
#########################

gridss_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/gridss",pattern = "*.bedpe", full.names= TRUE) 

gridss<-lapply(gridss_file, function(x){

    bed<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(gridss)){
    gridss[[i]]<-cbind(gridss[[i]],gridss_file[i])
    }

gridss_illu<-do.call("rbind",gridss)
names(gridss_illu)[length(names(gridss_illu))]<-"path"

gridss_illu_good<-gridss_illu%>% 
    mutate(sample_caller=sub('.*/\\s*', '', gsub("_gridss_survivor.bedpe","",path)), )  %>% 
    dplyr::select(!(path)) %>% 
    #mutate(sample=gsub("_.*$","",sample_caller)) %>% 
    dplyr::mutate(across(V2:V3,~.-1)) %>% 
    dplyr::mutate(across(V5:V6,~.-1)) 
    #mutate(unique_id =paste(sample,V11, sep = "_"))  ### V11 is the SV type
   
    
colnames(gridss_illu_good) <- c("chromA","startA","endA","chromB","startB","endB","typeID","coma","strand1","strand2","SVtype","sample") 
#chroms_good<-c(paste("chr",seq(1:22),sep=""),"chrX","chrY")

gridss_illu_good_ch<-gridss_illu_good[(gridss_illu_good$chromA %in% pan_chroms) & (gridss_illu_good$chromB %in% pan_chroms),]  ### remove unplaced scaffolds


### find overlaps     
samples<-unique(pan_uniq$sample)  ### what are the samples in the panel 



for(ss in 1:length(samples)){
  
    ## SV callers from resequenced 
    focal_caller<-gridss_illu_good_ch %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(chromA,startA,startB,SVtype,sample)
    #setDT(focal_caller_format)
    #setkey(focal_caller_format, chrom1,start1,start2)
    
    ## longread panel
    focal_panel<-pan_uniq %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(ChrA,ChrPosA,ChrPosB,SVtype) 

    colnames(focal_panel)<-c("chromA","startA","startB","SVtype")
    
    #rename_with(pan_good_ch, recode, ChrA = "chrom1",ChrPosA ="start1" ,  ChrPosB="start2" , SVtype= "SVtype")  ### fix renames


    chroms<-unique(focal_panel$chromA)
    
    out_res_chrom<-NULL
    for(ii in 1:length(chroms)){
        #for(ii in 1:20){

  
        ## SV caller
        focal_caller_chrom<-focal_caller %>%
         #filter(chromA ==chroms[ii] & SVtype != "TRA") %>%  ### tTO FIX: emporarily skip TRA
         filter(chromA ==chroms[ii]) %>% 
         filter(!(startB < startA))              ### tTO FIX: emporarily skip TRA
         
           
        
        ## panel
        focal_panel_chrom<-focal_panel %>% 
         filter(chromA ==chroms[ii])%>%  ### tTO FIX: emporarily skip TRA
         #filter(chromA ==chroms[ii] & SVtype != "TRA")%>%  ### tTO FIX: emporarily skip TRA
         filter(!(startB < startA))

        setDT(focal_caller_chrom)
        setDT(focal_panel_chrom)

        
        
        setkey(focal_caller_chrom, chromA,startA,startB)
        setkey(focal_panel_chrom, chromA,startA,startB)

        overlaps<-data.frame(foverlaps(focal_caller_chrom, focal_panel_chrom, type="any"))
        #overlaps_matchSV<-overlaps[overlaps$SVtype == overlaps$i.SVtype,] %>% 
        overlaps_matchSV<-overlaps %>% 
           drop_na()

        perfect_match<-overlaps_matchSV[overlaps_matchSV$startA==overlaps_matchSV$i.startA & overlaps_matchSV$startB==overlaps_matchSV$i.startB,]   


           ### rbind chroms
       if (nrow(perfect_match)>0) {

            out_res_chrom<-rbind(perfect_match,out_res_chrom) 


       } ## if loop 


    }  ### chrom loop 

}    

out_res_chrom<-out_res_chrom[,c(1,2,3,4,5,6,8)]
write.table(out_res_chrom, file ="~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/outputs/out_res_longread_pannel_resequenced_short_read_SW_Gridss.table", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


#########################
##################
### from here ###
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
