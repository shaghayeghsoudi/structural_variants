

### find intersection between selected panel from longread and short read SV ##
rm(list = ls())

### load required libraries
library(dplyr)
#library(ggplot2)
#library(vcfR)
#library(StructuralVariantAnnotation)
#library(plyr)
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

pan_chroms<-unique(c(pan_good$ChrA ,pan_good$ChrB))  ### unique chromosmes present in the PacBio pannel


###################################################
### load bedpe files (detected by just gridss) ###
###################################################

gridss_file<-list.files("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/gridss/gridss_annotated_bed_1index",pattern = "*.bed", full.names= TRUE) 

gridss<-lapply(gridss_file, function(x){

    bed<-read.delim(x, header = FALSE, sep = "\t")
})

for (i in 1:length(gridss)){
    gridss[[i]]<-cbind(gridss[[i]],gridss_file[i])
    }

gridss_illu<-do.call("rbind",gridss)
names(gridss_illu)[length(names(gridss_illu))]<-"path"

gridss_illu_good<-gridss_illu%>% 
    mutate(sample_caller=sub('.*/\\s*', '', gsub(".annotated.simple.bed","",path)), )  %>% 
    dplyr::select(!(path)) %>% 
    mutate(sample=gsub("_.*$","",sample_caller)) %>% 
    dplyr::select(c(V1:V5,sample))
    #mutate(sample=gsub("_.*$","",sample_caller)) %>% 
    #dplyr::mutate(across(V2:V3,~.-1),across(V5:V6,~.-1)) 
    #mutate(unique_id =paste(sample,V11, sep = "_"))  ### V11 is the SV type
colnames(gridss_illu_good) <- c("chromA","startA","endA","SVType","QUAL", "sample") 


### duplicates columns to match for datatable
gridss_illu_good_duplicated<-gridss_illu_good %>% 
    mutate(chromB = chromA, startB= startA, endB=endA) %>% 
    mutate(uniq_identifier=paste(sample,SVType,chromA,startA,chromB,endB, sep = "_"))%>% 
    mutate(sample_caller = paste(sample,"gridss",sep ="_"))
gridss_illu_good_duplicated <- gridss_illu_good_duplicated[,c("chromA","startA","endA","chromB","startB","endB","SVType","QUAL", "sample","uniq_identifier","sample_caller")]

gridss_illu_good_ch<-subset(gridss_illu_good_duplicated,chromA %in% pan_chroms & chromB %in% pan_chroms)
   
################################################################################################################
### load bedpe files (variants detected by each short read SV caller on re-sequenced cell lines/patient data) ###
################################################################################################################

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
    dplyr::mutate(across(V2:V3,~.-1),across(V5:V6,~.-1)) %>% 
    mutate(uniq_identifier=paste(sample,V11,V1,V2,V4,V5, sep = "_"))
    #mutate(unique_id =paste(sample,V11, sep = "_"))  ### V11 is the SV type
   



colnames(beds_illu_good) <- c("chromA","startA","endA","chromB","startB","endB","typeID","strand1","strand2","SVtype","sample_caller","sample","uniq_identifier") 

#chroms_good<-c(paste("chr",seq(1:22),sep=""),"chrX","chrY")

#beds_illu_good_ch<-beds_illu_good[(beds_illu_good$chromA %in% pan_chroms) & (beds_illu_good$chromB %in% pan_chroms),]  ### remove unplaced scaffolds
beds_illu_good_ch<-subset(beds_illu_good,chromA %in% pan_chroms & chromB %in% pan_chroms)

##########################################################################
### find overlaps  between long read panel and short read resequenced ####
##########################################################################

samples<-unique(pan_good$sample)  ### what are the samples in the panel 


for(ss in 1:length(samples)){
  
    ## SV callers from resequenced 
    focal_caller<-beds_illu_good_ch %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(chromA,startA,startB,SVtype,sample_caller,uniq_identifier)
    #setDT(focal_caller_format)
    #setkey(focal_caller_format, chrom1,start1,start2)
    
    ## longread panel
    focal_panel<-pan_good %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(ChrA,ChrPosA,ChrPosB,SVtype,uniq_identifier) 
    colnames(focal_panel)<-c("chromA","startA","startB","SVtype","uniq_identifier")
    
    #rename_with(pan_good_ch, recode, ChrA = "chrom1",ChrPosA ="start1" ,  ChrPosB="start2" , SVtype= "SVtype")  ### fix renames
    
    chroms<-unique(focal_panel$chromA)
    
    out_res_chrom<-NULL
    for(ii in 1:length(chroms)){
        #for(ii in 1:18){

        ## SV caller
        focal_caller_chrom<-focal_caller %>%
        filter(chromA ==chroms[ii] & SVtype != "TRA") %>%   ### toDO FIX: emporarily skip TRA
        filter(!(startB < startA))    
        setDT(focal_caller_chrom) 
        

         focal_caller_chrom_tra<-  focal_caller %>%  ### just TRA
        filter(chromA ==chroms[ii] & SVtype == "TRA")          
         
        ## panel
        focal_panel_chrom<-focal_panel %>% 
        filter(chromA ==chroms[ii] & SVtype != "TRA")%>%  ### tTO FIX: emporarily skip TRA
        filter(!(startB < startA))
        setDT(focal_panel_chrom)


        focal_panel_tra<-  focal_panel %>%  ### just TRA
        filter(chromA ==chroms[ii] & SVtype == "TRA")  

        

        setkey(focal_caller_chrom, chromA,startA,startB)
        setkey(focal_panel_chrom, chromA,startA,startB)

        overlaps<-data.frame(foverlaps(focal_caller_chrom, focal_panel_chrom, type="any"))
        overlaps_full_partial_matchSV<-overlaps[overlaps$SVtype == overlaps$i.SVtype,] %>%   ### match in SV type between panel and resequenced variants
           drop_na() %>% 
           mutate(start_diff = abs(startA-i.startA), end_diff = abs(startB-i.startB)) %>% 
           mutate(short_read_caller_chrom=chromA) %>% 
           dplyr::select(chromA,startA,startB,SVtype,short_read_caller_chrom,i.startA,i.startB,i.SVtype,sample_caller,start_diff,end_diff) 


            #colnames(overlaps_full_partial_matchSV)<-c("chrom","panel_start","panel_end","SV_type","short_read_caller_start","short_read_caller_end","sample_caller", "start_diff", "end_diff")
        others<-overlaps_full_partial_matchSV %>% 
        filter(!(SVtype == "INS") & SVtype==i.SVtype &start_diff <=50 & end_diff <=50)

           
       just_ins<-overlaps_full_partial_matchSV %>% 
       filter(SVtype=="INS" & start_diff <=50 & SVtype==i.SVtype )

       all_types<-rbind(others,just_ins) ### INS, DEL, DUP, INV
       all_types<-all_types[,c("chromA",   "startA" , "startB" ,"SVtype"  ,   "short_read_caller_chrom", "i.startA" ,"i.startB", "i.SVtype" ,"sample_caller" ,"start_diff", "end_diff")]
       colnames(all_types)<-c("panel_chrom","panel_start","panel_end","panel_SV_type","short_read_caller_chrom","short_read_caller_start","short_read_caller_end","short_red_caller_SV_type","sample_caller" ,"start_diff", "end_diff")
        
        ####################################
        ### join TRA panel and illumina ####
        ####################################
        join_tra<-merge(focal_panel_tra,focal_caller_chrom_tra,by.x = "uniq_identifier", by.y ="uniq_identifier" )
        join_tra
        if (nrow(join_tra > 0)){

           join_tra<-join_tra[,c("chromA.x","startA.x", "startB.x", "SVtype.x" ,"chromA.y" , "startA.y","startB.y", "SVtype.y","sample_caller")]
           colnames(join_tra)<-c("panel_chrom","panel_start","panel_end","panel_SV_type","short_read_caller_chrom","short_read_caller_start","short_read_caller_end","short_red_caller_SV_type","sample_caller")
           join_tra$start_diff<-0
           join_tra$end_diff<-0
        }

         tra_overlaps_full_partial_matchSV<-rbind(all_types,join_tra)
        

        #perfect_match<-overlaps_matchSV[overlaps_matchSV$startA==overlaps_matchSV$i.startA & overlaps_matchSV$startB==overlaps_matchSV$i.startB,]   


           ### rbind chroms
       if (nrow(tra_overlaps_full_partial_matchSV)>0) {

            out_res_chrom<-rbind(tra_overlaps_full_partial_matchSV,out_res_chrom) 


       } ## if loop 


    }  ### chrom loop 
    write.table(out_res_chrom, file =paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/outputs-celllines_pacbiopanel_ovelap_resequenced_illumina_full_partial_overlap/out_res_PacBio_pannel_overlap_resequenced_illumina_",samples[ss],"_other3_full_partial_50base.table", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


}    


################################################################
###### ### find overlaps between GRIDSS and PacBio panel #######
################################################################



for(ss in 1:length(samples)){
  
    ## SV callers from resequenced 
    focal_girdss<-gridss_illu_good_ch %>% 
    filter(sample==samples[ss]) %>% 
    #dplyr::select(chromA,startA,startB,SVtype,sample)
    dplyr::select(chromA,startA,endA,SVType,sample) %>% 
    mutate(sample = paste(sample,"gridss", sep = "_"))

    #setDT(focal_caller_format)
    #setkey(focal_caller_format, chrom1,start1,start2)
    
    ## longread panel
    focal_panel<-pan_good %>% 
    filter(sample==samples[ss]) %>% 
    dplyr::select(ChrA,ChrPosA,ChrPosB,SVtype) 

    colnames(focal_panel)<-c("chromA","startA","startB","SVType")
    
    #rename_with(pan_good_ch, recode, ChrA = "chrom1",ChrPosA ="start1" ,  ChrPosB="start2" , SVtype= "SVtype")  ### fix renames


    chroms<-unique(focal_panel$chromA)
    
    out_res_chrom_gridss<-NULL
    for(ii in 1:length(chroms)){
        #for(ii in 1:20){

         focal_gridss_chrom<-focal_girdss %>%
         #filter(chromA ==chroms[ii] & SVtype != "TRA") %>%  ### tTO FIX: emporarily skip TRA
         filter(chromA ==chroms[ii]) %>% 
         filter(!(startA > endA))              ### tTO FIX: emporarily skip TRA
         
        
        ## PacBio panel
        focal_panel_chrom<-focal_panel %>% 
         filter(chromA ==chroms[ii])%>%  ### tTO FIX: emporarily skip TRA
         #filter(chromA ==chroms[ii] & SVtype != "TRA")%>%  ### tTO FIX: emporarily skip TRA
         filter(!(startB < startA))

        setDT(focal_gridss_chrom)
        setDT(focal_panel_chrom)

        setkey(focal_gridss_chrom, chromA,startA,endA)
        setkey(focal_panel_chrom, chromA,startA,startB)

        overlaps_gridss<-data.frame(foverlaps(focal_gridss_chrom, focal_panel_chrom, type="any"))
        #overlaps_gridss<-overlaps_gridss[overlaps_gridss$SVtype == overlaps_gridss$i.SVtype,]


        overlaps_matchSV_gridss<-overlaps_gridss %>% 
        drop_na() %>% 
        mutate(start_diff = abs(startA-i.startA), end_diff = abs(startB-endA)) 

        
        others<-overlaps_matchSV_gridss %>% 
        filter(!(SVType == "INS") & SVType==i.SVType &start_diff <=50 & end_diff <=50)

           
       just_ins<-overlaps_matchSV_gridss %>% 
       filter(SVType=="INS" & start_diff <=50 & SVType==i.SVType )

       all_types<-rbind(others,just_ins)

        #colnames(overlaps_matchSV_gridss)<-c("chrom","panel_start","panel_end","SV_type","short_read_caller_start","short_read_caller_end","sample_caller", "start_diff", "end_diff")


           ### rbind chroms
       if (nrow( all_types)>0) {

             out_res_chrom_gridss<-rbind(all_types, out_res_chrom_gridss) 


       } ## if loop 


    }  ### chrom loop 

    write.table(out_res_chrom_gridss, file =paste("~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/outputs-celllines_pacbiopanel_ovelap_resequenced_illumina_full_partial_overlap/out_res_PacBio_pannel_overlap_resequenced_illumina_",samples[ss],"_gridss_full_partial_50base.table", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE)


}    


#########################
###############
### create matrix
#results_out <- array (NA, c((nrow (input_good) * (end1 - start1 + 1)),7))
#count <- 0

#system.time(
##loop through focal snp
#for (i in start1:end1){
#	
#	loop through all SNPs
#	for (j in 1:nrow (input_good)){
#		count <- count + 1
#		 results_out[count,7] <- cor (input_good[i,], input_good[j,], use = "pairwise.complete.obs")
#		results_out [count,1] <- scafpos_good[i,1]
#		results_out [count,2] <- scafpos_good[i,2]
#		#results_out [count,3] <- scafpos_good[i,3]
#		results_out [count,4] <- scafpos_good[j,1]
#		results_out [count,5] <- scafpos_good[j,2]
#		#results_out [count,6] <- scafpos_good[j,3]
#	
#	}
#}
#)


#outname <- paste ("output", start1,"_",end1,".txt",sep = "")
#write.table (results_out, outname, col.names = F, row.names = F, quote = F)
