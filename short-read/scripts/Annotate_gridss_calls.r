### Annotate gridss calls with SIMPLE_TYPE and SVLEN fields

 ## NOTE: Please adjust the name of your vcf file  in the script!!!

#BiocManager::install("StructuralVariantAnnotation")
#install.packages("stringr")
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
#' Simple SV type classifier
simpleEventType <- function(gr) {
  pgr = partner(gr)
  return(ifelse(seqnames(gr) != seqnames(pgr), "CTX", # inter-chromosomosal
    ifelse(strand(gr) == strand(pgr), "INV",
      ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
        ifelse(xor(start(gr) < start(pgr), strand(gr) == "-"), "DEL",
          "DUP")))))
}
# using the example in the GRIDSS /example directory
root<-"~/Dropbox/cancer_reserach/sarcoma/sarcoma_analysis/short_reads_SV/cell_lines_resequenced/pacbio_resequenced_short_read/gridss/"
vcf <- readVcf(paste(root,"GCT_gridss.filter.pass.vcf", sep = ""), "hg19")         ##### NOTE: change here based on the name of your vcf file
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
	row.names=c("SIMPLE_TYPE"),
	Number=c("1"),
	Type=c("String"),
	Description=c("Simple event type annotation based purely on breakend position and orientation."))), "DataFrame"))
gr <- breakpointRanges(vcf)

### error message pops up: Warning message:
#In .breakpointRanges(x, ...) : Removing 12 unpaired breakend variants.
## Use inferMissingBreakends=TRUE to recover with inferred partner breakends.
#gr <- breakpointRanges(vcf,inferMissingBreakends=TRUE)
 
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf[gr$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[gr$sourceId])$SVLEN <- gr$svLen

 ##### NOTE: change here based on the name of your vcf file
writeVcf(vcf, "chr12.1527326.DEL1024.sv.annotated.vcf") # generated by example/gridss.sh  

# TODO: perform event filtering here
# By default, GRIDSS is very sensitive but this comes at the cost of a high false discovery rate
gr <- gr[gr$FILTER == "PASS" & partner(gr)$FILTER == "PASS"] # Remove low confidence calls

simplegr <- gr[simpleEventType(gr) %in% c("INS", "INV", "DEL", "DUP")]
simplebed <- data.frame(
	chrom=seqnames(simplegr),
	# call the centre of the homology/inexact interval
	start=as.integer((start(simplegr) + end(simplegr)) / 2),
	end=as.integer((start(partner(simplegr)) + end(partner(simplegr))) / 2),
	name=simpleEventType(simplegr),
	score=simplegr$QUAL,
	strand="."
)
# Just the lower of the two breakends so we don't output everything twice
simplebed <- simplebed[simplebed$start < simplebed$end,]


### change here based on the name of vcf file
write.table(simplebed, "chr12.1527326.DEL1024.simple.bed", quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)