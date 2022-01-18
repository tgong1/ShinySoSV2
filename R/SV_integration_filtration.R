#' Convert VCF format to bed format
#'
#' This function read vcf file and convert it to bed format
#'
#' @param vcf_file names of vcf file
#' @return data frame
#' @export
vcf_to_bed <- function(vcf_file){
  vcf <- VariantAnnotation::readVcf(vcf_file)

  fixed_df <- vcf@fixed
  gr <- vcf@rowRanges
  info <- vcf@info

  for(field in c("SVLEN","END", "STRANDS","CT","INV5","INV3","MATEID")){
    tmp <- eval(parse(text=paste0("info$", field)))
    if (length(tmp) == 0){
      assign(paste0("INFO_", field), NA)
    }else{
      idx <- !(sapply(tmp, length))
      tmp[idx] <- NA
      if(field == "STRANDS"){
        assign(paste0("INFO_", field), sapply(tmp, paste, collapse = ","))
      }else{
        assign(paste0("INFO_", field), unlist(tmp))
      }
    }
  }

  bed <- data.frame(CHROM = gr@seqnames,
                    POS = gr@ranges@start,
                    ID_caller = gr@ranges@NAMES,
                    REF = as.character(fixed_df$REF),
                    ALT = data.frame(fixed_df$ALT)$value,
                    QUAL = fixed_df$QUAL,
                    FILTER = fixed_df$FILTER,
                    INFO_END,
                    INFO_SVTYPE = info$SVTYPE,
                    INFO_SVLEN,
                    INFO_STRANDS,
                    INFO_CT = INFO_CT,
                    INFO_INV5 = INFO_INV5,
                    INFO_INV3 = INFO_INV3,
                    INFO_MATEID_caller = INFO_MATEID,
                    stringsAsFactors = FALSE
  )

  return(bed)
}

#' Classify SV types based on VCF
#'
#' This function read bed format
#'
#' @param bed datafame of bed format
#' @param caller_name name of caller to be used in ID
#' @return data frame
#' @export
simple_SVTYPE_classification <- function(bed, caller_name){
  bedpe <- bed_to_bedpe(bed)
  if(length(bedpe$ID_caller)!=0){
    bedpe <- bedpe[is.na(match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)) | ### either don't have mate (i.e. not BND)
                     (c(1:nrow(bedpe)) < match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)),] ###  or present first as BND
  }

  SV_index <- c(1:nrow(bedpe))
  event_index <- SV_index

  SVTYPE <- bedpe$INFO_SVTYPE
  SVTYPE[(bedpe$strand1 == "+") & (bedpe$strand2 == "-") & (bedpe$chrom1 == bedpe$chrom2)] <- "DEL"
  SVTYPE[(bedpe$strand1 == "-") & (bedpe$strand2 == "+") & (bedpe$chrom1 == bedpe$chrom2)] <- "DUP"

  SVTYPE[(bedpe$strand1 == bedpe$strand2) & (bedpe$chrom1 == bedpe$chrom2)] <- "INV"
  SVTYPE[(bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  SVTYPE[is.na(bedpe$strand1)] <- "INS"

  standard_bedpe <- data.frame(chrom1 = as.character(bedpe$chrom1),
                               pos1 = as.integer(bedpe$pos1),
                               chrom2 = as.character(bedpe$chrom2),
                               pos2 = as.integer(bedpe$pos2),
                               SVTYPE = SVTYPE,
                               strand1 = bedpe$strand1,
                               strand2 = bedpe$strand2,
                               ID = paste0(caller_name,"_",SV_index,"_1_",event_index),
                               ID_mate = paste0(caller_name,"_",SV_index,"_2_",event_index),
                               #ALT = bedpe$ALT,
                               stringsAsFactors = FALSE)

  bedpe_SVTYPE_classified <- data.frame(standard_bedpe, bedpe[,!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2"))])
  colnames(bedpe_SVTYPE_classified) <- c("chrom1", "pos1","chrom2","pos2","SVTYPE",
                                         "strand1","strand2",
                                         "ID","ID_mate",
                                         colnames(bedpe)[!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2"))])
  return(bedpe_SVTYPE_classified)
}
