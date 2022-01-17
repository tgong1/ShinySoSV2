#' Convert bed format to bedpe format
#'
#' This function read bed format
#'
#' @param bed dataframe
#' @return data frame
#' @export
bed_to_bedpe <- function(bed){
  if(sum(!is.na(bed$INFO_STRANDS))!=0){
    tmp <- data.frame(str_split_fixed(bed$INFO_STRANDS,",",2), bed$ID_caller)
    tmp$X1 <- gsub(":.*","", tmp$X1)
    tmp$X2 <- gsub(":.*","", tmp$X2)
    idx <- rep(1:nrow(tmp), rowSums(tmp[,c(1,2)]!=""))
    tmp2 <- tmp[idx,]

    tmp2[duplicated(tmp2$bed.ID_caller),]$X1 <- tmp2[duplicated(tmp2$bed.ID_caller),]$X2
    tmp3 <- tmp2[,-2]
    tmp3 <- cbind(tmp3,str_split_fixed(tmp3$X1,"",2))
    colnames(tmp3) <- c("strands","ID_caller","strand1","strand2")
    strand1 <- tmp3$strand1
    strand2 <- tmp3$strand2
    #bed <- bed[match(tmp3$ID_caller, bed$ID),]
    #bed$strand1 <- tmp3$strand1
    #bed$strand2 <- tmp3$strand2
  }
  if(sum(!is.na(bed$INFO_CT)) != 0){
    strand1 <- rep(NA, nrow(bed)); strand2 <- rep(NA, nrow(bed))
    strand1[bed$INFO_CT=="3to5"] <- "+"; strand2[bed$INFO_CT=="3to5"] <- "-"
    strand1[bed$INFO_CT=="5to3"] <- "-"; strand2[bed$INFO_CT=="5to3"] <- "+"
    strand1[bed$INFO_CT=="3to3"] <- "+"; strand2[bed$INFO_CT=="3to3"] <- "+"
    strand1[bed$INFO_CT=="5to5"] <- "-"; strand2[bed$INFO_CT=="5to5"] <- "-"
    strand1[bed$INFO_CT=="NtoN"] <- NA; strand2[bed$INFO_CT=="NtoN"] <- NA
  }

  if(sum(!is.na(bed$INFO_INV5)) !=0){
    strand1 <- rep(NA, nrow(bed)); strand2 <- rep(NA, nrow(bed))
    strand1[bed$INFO_INV5] <- "-";strand2[bed$INFO_INV5] <- "-"
    strand1[bed$INFO_INV3] <- "+";strand2[bed$INFO_INV3] <- "+"
  }

  ALT <- bed$ALT
  tmp <- gsub("\\:.*",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub(".*\\[",'', tmp)
  tmp2 <- gsub(".*\\]",'', tmp1)
  chrom1 = as.character(bed$CHROM)
  chrom2 <- chrom1
  chrom2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2

  tmp <- gsub(".*\\:",'', ALT[grepl('\\[',ALT) | grepl('\\]',ALT)])
  tmp1 <- gsub("\\[.*",'', tmp)
  tmp2 <- gsub("\\].*",'', tmp1)

  pos2 <- bed$INFO_END
  pos2[grepl('\\[',ALT) | grepl('\\]',ALT)] <- tmp2
  pos2[is.na(pos2)] <- bed$POS[is.na(pos2)] ###if no POS2 then use #CHROM in VCF



  bedpe <- data.frame(chrom1=bed$CHROM,
                      pos1=bed$POS,
                      chrom2=paste0(tolower(substring(chrom2,1,3)), substring(chrom2,4,nchar(chrom2))),
                      pos2=pos2,
                      ALT=ALT,
                      stringsAsFactors = FALSE)

  bedpe <- data.frame(bedpe, bed[,!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT"))])
  colnames(bedpe) <- c("chrom1", "pos1","chrom2","pos2","ALT", colnames(bed)[!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT"))])
  return(bedpe)
}
