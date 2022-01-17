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
  if(!isEmpty(bedpe$ID)){
    bedpe <- bedpe[is.na(match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)) | ### either don't have mate (i.e. not BND)
                     (c(1:nrow(bedpe)) < match(bedpe$ID_caller, bedpe$INFO_MATEID_caller)),] ###  or present first as BND
  }

  SV_index <- c(1:nrow(bedpe))
  #ID_tmp <- bedpe$EVENT
  #ID_tmp[is.na(bedpe$EVENT)] <- SV_index[is.na(bedpe$EVENT)]
  #event_index <- match(ID_tmp,unique(ID_tmp))
  event_index <- SV_index

  strand1 <- rep(NA,nrow(bedpe))
  strand2 <- rep(NA,nrow(bedpe))

  strand1[bedpe$INFO_SVTYPE == "DEL"] <- "+"
  strand2[bedpe$INFO_SVTYPE == "DEL"] <- "-"
  strand1[bedpe$INFO_SVTYPE == "DUP"] <- "-"
  strand2[bedpe$INFO_SVTYPE == "DUP"] <- "+"

  ###[p[t
  strand1[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) == "["] <- "-"
  strand2[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) == "["] <- "-"

  ###t[p[
  strand1[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) != "["] <- "+"
  strand2[grepl('\\[', bedpe$ALT) & substr(bedpe$ALT,1,1) != "["] <- "-"

  ###t]p]
  strand1[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) != "]"] <- "+"
  strand2[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) != "]"] <- "+"

  ###]p]t
  strand1[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) == "]"] <- "-"
  strand2[grepl(']', bedpe$ALT) & substr(bedpe$ALT,1,1) == "]"] <- "+"


  SVTYPE <- bedpe$INFO_SVTYPE
  SVTYPE[(strand1 == "+") & (strand2 == "-") & (bedpe$chrom1 == bedpe$chrom2)] <- "DEL"
  SVTYPE[(strand1 == "-") & (strand2 == "+") & (bedpe$chrom1 == bedpe$chrom2)] <- "DUP"

  SVTYPE[(strand1 == strand2) & (bedpe$chrom1 == bedpe$chrom2)] <- "INV"
  #SVTYPE[(strand1 == strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA_INV"
  #SVTYPE[(strand1 != strand2) & (bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  SVTYPE[(bedpe$chrom1 != bedpe$chrom2)] <- "TRA"
  SVTYPE[is.na(strand1)] <- "INS"

  standard_bedpe <- data.frame(chrom1 = as.character(bedpe$chrom1),
                               pos1 = as.integer(bedpe$pos1),
                               chrom2 = as.character(bedpe$chrom2),
                               pos2 = as.integer(bedpe$pos2),
                               SVTYPE = SVTYPE,
                               #strand1 = strand1,
                               #strand2 = strand2,
                               ID = paste0(caller_name,"_",SV_index,"_1_",event_index),
                               ID_mate = paste0(caller_name,"_",SV_index,"_2_",event_index),
                               ALT = bedpe$ALT,
                               stringsAsFactors = FALSE)

  bedpe_SVTYPE_classified <- data.frame(standard_bedpe, bedpe[,!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2","ALT"))])
  colnames(bedpe_SVTYPE_classified) <- c("chrom1", "pos1","chrom2","pos2","SVTYPE",
                                         # "strand1","strand2",
                                         "ID","ID_mate","ALT",
                                         colnames(bedpe)[!(colnames(bedpe) %in% c("chrom1","chrom2","pos1","pos2","ALT"))])
  return(bedpe_SVTYPE_classified)
}

