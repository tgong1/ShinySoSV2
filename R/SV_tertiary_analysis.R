#' SVtype count
#'
#' This function read bed format
#'
#' @param bedpe bed format
#' @return data frame
#' @export
SVTYPE_stat_generate <- function(bedpe){
  All_SVTYPE <- unique(bedpe$SVTYPE)
  for(i in c(1: length(All_SVTYPE))){
    assign(paste0("N_", All_SVTYPE[i]), sum(bedpe$SVTYPE == All_SVTYPE[i]))
  }
  STAT_bed <- data.frame(do.call("cbind", lapply(paste0("N_", All_SVTYPE),function(s) eval(parse(text=s)))))
  colnames(STAT_bed) <- All_SVTYPE
  return(STAT_bed)
}

#' Summary SV count
#'
#' This function read bed format
#'
#' @param All_sampleID sample ID
#' @param All_input_df_name names
#' @return data frame
#' @export
Summary_SV_type <- function(All_sampleID, All_input_df_name){
  summary_results <- c()
  for(i in c(1:length(All_sampleID))){
    sampleID <- All_sampleID[i]
    df <- eval(parse(text = All_input_df_name[i]))
    SVTYPE_count <- SVTYPE_stat_generate(df)

    all_colnames <- unique(c(colnames(summary_results), colnames(SVTYPE_count)))

    SVTYPE_count_tmp <- data.frame(matrix(0, nrow = nrow(SVTYPE_count), ncol = length(all_colnames)))
    colnames(SVTYPE_count_tmp) <- all_colnames
    SVTYPE_count_tmp[colnames(SVTYPE_count_tmp) %in% colnames(SVTYPE_count)] <- SVTYPE_count

    if(i > 1){
      summary_results_tmp <- data.frame(matrix(0, nrow = nrow(summary_results), ncol = length(all_colnames)))
      colnames(summary_results_tmp) <- all_colnames
      summary_results_tmp[colnames(summary_results_tmp) %in% colnames(summary_results)] <- summary_results
      summary_results <- summary_results_tmp
    }
    summary_results <- rbind(summary_results, SVTYPE_count_tmp)
  }
  summary_results <- data.frame(sampleID =  All_sampleID, summary_results)
  return(summary_results)
}

#' Spectrum of SV count
#'
#' This function read bed format
#'
#' @param input_SV_count SV type count data frame
#' @param threshold_total threshold of minimum total count of SVs per sample
#' @param threshold_relative_freq threshold of minimum relative frequency of one SV type
#' @return data frame of hyper SV
#' @export
Spectrum_SV_type <- function(input_SV_count, threshold_total, threshold_relative_freq){
  theme1 <-  ggplot2::theme(axis.text=ggplot2::element_text(size=12,face="bold"),
                            axis.title=ggplot2::element_text(size=14,face="bold"),
                            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8),
                            plot.title = ggplot2::element_text(size=14),
                            legend.text = ggplot2::element_text(size=12,face="bold"),
                            #legend.title = element_text(size=12,face="bold"),
                            legend.title = ggplot2::element_blank(),
                            legend.position="top")

  theme2 <-  ggplot2::theme(axis.text = ggplot2::element_text(size=12,face="bold"),
                            axis.title=ggplot2::element_text(size=14,face="bold"),
                            plot.title = ggplot2::element_text(size=14),
                            legend.text = ggplot2::element_text(size=12,face="bold"),
                            #legend.title = element_text(size=12,face="bold"),
                            legend.title = ggplot2::element_blank(),
                            legend.position="top")
  df <- data.frame(SVTYPE = rep(colnames(input_SV_count)[2:ncol(input_SV_count)], each = nrow(input_SV_count)),
                   sampleID = rep(input_SV_count$sampleID, (ncol(input_SV_count)-1)),
                   Count = as.vector(unlist(input_SV_count[,2:ncol(input_SV_count)])))

  p1 <- ggplot2::ggplot(data=df, ggplot2::aes(x=sampleID, y=Count, fill=SVTYPE)) +
    ggplot2::geom_bar(stat="identity")+
    ggplot2::scale_fill_manual(values = RColorBrewer::brewer.pal(12, "Set3"), drop=F)+
    theme1
  pdf(file="./Spectrum_SV_across_sample.pdf", width=19, height=9)
  print(p1)
  dev.off()

  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = SVTYPE, y=Count)) +
    ggplot2::geom_boxplot()+
    theme2

  N_total <- rowSums(input_SV_count[,2:ncol(input_SV_count)])
  if(missing(threshold_total)){
    threshold_total <- mean(N_total)
  }
  if(missing(threshold_relative_freq)){
    threshold_relative_freq <- 0.5
  }

  df2 <- cbind(df,
               relative_freq = df$Count/(rep(N_total,(ncol(input_SV_count)-1))),
               N_total = rep(N_total,(ncol(input_SV_count)-1)))
  p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = SVTYPE, y=relative_freq)) +
    ggplot2::geom_boxplot()+
    ggplot2::ylab("Relative frequency")+
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1))+
    theme2

  figure1 <- ggpubr::ggarrange(p1, p2,
                               labels = c("A", "B"), heights = c(1,1),
                               ncol = 2, nrow = 1)

  pdf(file="./Spectrum_SV_across_SVTYPE.pdf",width=9,height=6)
  print(figure1)
  dev.off()

  hyper_SV <- df2[df2$relative_freq > threshold_relative_freq & df2$N_total > threshold_total,]
  return(hyper_SV)
}
