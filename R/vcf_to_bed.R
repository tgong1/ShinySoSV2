#' Convert VCF format to bed format
#'
#' This function read
#'
#' @param vcf_file names of vcf file
#' @return data frame
#' @export
vcf_to_bed <- function(vcf_file){
  vcf <- VariantAnnotation::readVcf(vcf_file)
  info = VariantAnnotation::info(vcf)
  gr <- rowRanges(vcf)

  for(field in c("SVTYPE", "SVLEN", "END","MATEID","STRANDS","CT","INV5","INV3")){
    tmp <- eval(parse(text=paste0("info$", field)))
    idx <- !(sapply(tmp, length))
    if (isEmpty(idx)){
      assign(paste0("INFO_", field), NA)
    }else{
      tmp[idx] <- NA
      if(field == "STRANDS"){
        assign(paste0("INFO_", field), sapply(tmp, paste, collapse = ","))
      }else{
        assign(paste0("INFO_", field), unlist(tmp))
      }
    }
  }

  bed <- data.frame(CHROM = seqnames(gr),
                    POS = start(gr),
                    INFO_END,
                    INFO_SVTYPE,
                    INFO_SVLEN,
                    INFO_STRANDS,
                    INFO_CT,
                    #INFO_INV5 = ifelse(INFO_INV5,"INV5",NA),
                    #INFO_INV3 = ifelse(INFO_INV3,"INV3",NA),
                    INFO_INV5 = INFO_INV5,
                    INFO_INV3 = INFO_INV3,
                    REF = gr$REF,
                    ALT = unlist(gr$ALT),

                    ID_caller = names(gr),
                    QUAL = gr$QUAL,
                    FILTER = gr$FILTER,
                    INFO_MATEID_caller = INFO_MATEID,
                    # EVENT=info$EVENT,
                    stringsAsFactors = FALSE)

  return(bed)
}
