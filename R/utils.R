#' Convert strands into standard form
#'
#' This function read bed format
#'
#' @param bed dataframe
#' @return data frame
#' @export
strands_standardisation <- function(bed){
if(length(bed$ID_caller) == 0){bed$ID_caller = c(1:nrow(bed))}
strand1 <- rep(NA, nrow(bed)); strand2 <- rep(NA, nrow(bed))
if(sum(!is.na(bed$INFO_STRANDS))!=0){
  tmp <- data.frame(stringr::str_split_fixed(bed$INFO_STRANDS,",",2), bed$ID_caller)
  tmp$X1 <- gsub(":.*","", tmp$X1)
  tmp$X2 <- gsub(":.*","", tmp$X2)
  idx <- rep(1:nrow(tmp), rowSums(tmp[,c(1,2)]!=""))
  tmp2 <- tmp[idx,]

  tmp2[duplicated(tmp2$bed.ID_caller),]$X1 <- tmp2[duplicated(tmp2$bed.ID_caller),]$X2
  tmp3 <- tmp2[,-2]
  tmp3 <- cbind(tmp3, stringr::str_split_fixed(tmp3$X1,"",2))
  colnames(tmp3) <- c("strands","ID_caller","strand1","strand2")
}else if(sum(!is.na(bed$INFO_CT)) != 0){
  tmp <- data.frame(bed$INFO_CT, ID_caller = bed$ID_caller)
  strand1[tmp$INFO_CT=="3to5"] <- "+"; strand2[tmp$INFO_CT=="3to5"] <- "-"
  strand1[tmp$INFO_CT=="5to3"] <- "-"; strand2[tmp$INFO_CT=="5to3"] <- "+"
  strand1[tmp$INFO_CT=="3to3"] <- "+"; strand2[tmp$INFO_CT=="3to3"] <- "+"
  strand1[tmp$INFO_CT=="5to5"] <- "-"; strand2[tmp$INFO_CT=="5to5"] <- "-"
  strand1[tmp$INFO_CT=="NtoN"] <- NA; strand2[tmp$INFO_CT=="NtoN"] <- NA
  tmp3 <- cbind(tmp, strand1, strand2)
}else if(sum(!is.na(bed$INFO_INV5)) !=0){
  tmp <- data.frame(bed$INFO_INV5, bed$INFO_INV3, ID_caller = bed$ID_caller)
  strand1[bed$INFO_INV5] <- "-";strand2[bed$INFO_INV5] <- "-"
  strand1[bed$INFO_INV3] <- "+";strand2[bed$INFO_INV3] <- "+"
  tmp3 <- cbind(tmp, strand1, strand2)
}else{
  tmp3 <- cbind(bed, strand1, strand2)
}

bed2 <- bed[match(tmp3$ID_caller, bed$ID_caller),]
bed2$strand1 <- tmp3$strand1
bed2$strand2 <- tmp3$strand2

if((sum(bed2$INFO_SVTYPE %in% c("DEL","DUP")) != 0) & sum(is.na(bed2$strand1))!=0){
  bed2[is.na(bed2$strand1) & bed2$INFO_SVTYPE == "DEL",]$strand1 <- "+"
  bed2[is.na(bed2$strand2) & bed2$INFO_SVTYPE == "DEL",]$strand2 <- "-"
  bed2[is.na(bed2$strand1) & bed2$INFO_SVTYPE == "DUP",]$strand1 <- "-"
  bed2[is.na(bed2$strand2) & bed2$INFO_SVTYPE == "DUP",]$strand2 <- "+"
}
return(bed2)
}

#' Convert bed format to bedpe format
#'
#' This function read bed format
#'
#' @param bed data frame
#' @return data frame
#' @export
bed_to_bedpe <- function(bed){
  bed <- strands_standardisation(bed)

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

  strand1 <- bed$strand1; strand2 <- bed$strand2
  ###[p[t
  strand1[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) == "["] <- "-"
  strand2[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) == "["] <- "-"

  ###t[p[
  strand1[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) != "["] <- "+"
  strand2[grepl('\\[', bed$ALT) & substr(bed$ALT,1,1) != "["] <- "-"

  ###t]p]
  strand1[grepl(']', bed$ALT) & substr(bed$ALT,1,1) != "]"] <- "+"
  strand2[grepl(']', bed$ALT) & substr(bed$ALT,1,1) != "]"] <- "+"

  ###]p]t
  strand1[grepl(']', bed$ALT) & substr(bed$ALT,1,1) == "]"] <- "-"
  strand2[grepl(']', bed$ALT) & substr(bed$ALT,1,1) == "]"] <- "+"

  bedpe <- data.frame(chrom1=bed$CHROM,
                      pos1=bed$POS,
                      chrom2=paste0(tolower(substring(chrom2,1,3)), substring(chrom2,4,nchar(chrom2))),
                      pos2=pos2,
                      strand1 = strand1,
                      strand2 = strand2,
                      ALT=ALT,
                      stringsAsFactors = FALSE)

  bedpe <- data.frame(bedpe, bed[,!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT","strand1","strand2"))])
  colnames(bedpe) <- c("chrom1", "pos1","chrom2","pos2","strand1","strand2","ALT",
                       colnames(bed)[!(colnames(bed) %in% c("CHROM","POS","INFO_END","ALT","strand1","strand2"))])
  return(bedpe)
}
