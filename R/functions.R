
#' Title
#'
#' @param se
#'
#' @return
#' @export
#'
#' @examples
missing_values_visualization <- function (se) {

#Create a new data frame containing only peptides with missing values
df_missing <- SummarizedExperiment::assay({{se}}) %>%
    dplyr::filter(!complete.cases(.)) %>% #Removes rows (peptides) with no missing values
    dplyr::mutate_all(~ ifelse(is.na(.), 0, 1)) #NA's replaced with 0 and everything else with 1

#For annotation purposes, columns must be in a specific order.
col_order <- metadata %>%
    dplyr::arrange(time_num, site) %>%  #Arrange rows based on time and site
    row.names() #Keep a vector of the row names (which are our column names in df_missing)

#Now rearrange columns in df_missing
df_missing <- df_missing %>%
    dplyr::select(all_of(col_order))

#Heatmap
heatmap_missing <- pheatmap(df_missing,
                            cluster_rows = F,
                            cluster_cols = F,
                            annotation_col = dplyr::select(metadata, c("site", "time")),
                            annotation_colors=list(time=c(baseline="#003f5c",exercise="#58508d", pre="#bc5090", rec180="#ff6361", rec60="#ffa600"),
                                                   site=c(a="#f1010e", v="#799fcb")),
                            show_rownames = F,
                            color=colorRampPalette(c("white", "black"))(2),
                            #cellwidth =5,
                            border_color = NA,
                            legend = F,
                            annotation_legend = T
)
}

#' Check if any peptides match selected thresholds of missing values on artery and vein site
#'
#' @param se
#' @param timepoint
#' @param artery_missing_threshold
#' @param vein_missing_threshold
#'
#' @return Data table
#' @export
#'
#' @examples
check_missing_thresholds <- function (se, timepoint, artery, vein) {

    #Subset SE for only exercise values
    df <- SummarizedExperiment::assay({{se}}[, {{se}}$time=={{timepoint}}])
    se <- PhosphoExperiment(assay = list(abundance = df), colData=metadata[metadata$time == {{timepoint}},])

    #Substitute valid values with TRUE and missing with FALSE in both artery and vein matrix
    artery_presence <- !is.na(SummarizedExperiment::assay(se[,se$site == 'a']))

    vein_presence <- !is.na(SummarizedExperiment::assay(se[,se$site == 'v']))

    #Calculate percentages
    vein_detection_percentage <- rowMeans(vein_presence) * 100

    vein_detection_percentage <- vein_detection_percentage %>%
        as.data.frame()

    artery_detection_percentage <- rowMeans(artery_presence) * 100

    artery_detection_percentage <- artery_detection_percentage %>%
        as.data.frame()

    #Combine dataframes
    peptide_detection_percentage <- merge(vein_detection_percentage, artery_detection_percentage, by="row.names") %>%
        dplyr::rename(vein = "..x",
                      artery = "..y",
                      peptide = "Row.names")

    #Create data frame with peptides of interest (e.g., peptides with >70% valid in vein and 0% in artery)
    peptides_of_interest <- peptide_detection_percentage %>%
        dplyr::filter(artery <{{artery}} & vein > {{vein}})

    #Peptides that fit selected thresholds
    peptides_of_interest


}



## Below is a function to generate a volcano plot. It takes 3 arguments. "type" is the type of p-value adjustment
#(in the below code, xiao adjustment and q-values are added as columns, and thus arguments can be "xiao" or "q". "threshold" is the p-value or significance score cut-off.)

volcano <- function (dataset, type, threshold) {

    dataset%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Upregulated",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Downregulated",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=peptide))+
        geom_point(aes(color = color), size = 3)+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_blank(),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15, family="Source Sans Pro"),
              legend.title = element_blank(),
              legend.key = element_blank(),
              )+
        geom_text_repel(point.size=4, size=3, min.segment.length = Inf, force=0.3)+
        scale_color_manual(breaks = c("Upregulated", "Downregulated", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        xlab("Log2fold change") + ylab("-log10(p)")
}


volcano_all <-function (type, threshold) {

    all_results <- rbind(results_one_sample_pre, results_one_sample_baseline, results_one_sample_exercise, results_one_sample_rec60, results_one_sample_rec180)

    all_results%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Uptake",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Release",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=metabolite))+
        geom_point(aes(color = color), size = 3)+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_rect(color="white"),
              plot.title = element_text(size=12),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15), #family="Source Sans Pro"),
              legend.title = element_blank(),
              legend.key = element_blank(),
        )+
        geom_text_repel(point.size=4, size=3, min.segment.length = Inf, force=0.3)+
        scale_color_manual(breaks = c("Uptake", "Release", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        xlab("Log2fold change") + ylab("-log10(p)")+
        facet_grid(~factor(time, levels=c('pre', 'baseline', 'exercise', 'rec60', 'rec180')))

}

volcano_all_two_sample <-function (type, threshold) {

    all_results <- rbind(results_limma_pre, results_limma_exercise, results_limma_rec60, results_limma_rec180)

    all_results%>%
        mutate(color = case_when(
            logFC >= 0 & {{type}} <= {{threshold}} ~ "Uptake",
            logFC <= 0 & {{type}} <= {{threshold}} ~ "Release",
            TRUE ~ "Unchanged")) %>%
        ggplot(aes(x=logFC, y=-log10(P.Value), label=metabolite))+
        geom_point(aes(color = color), size = 3)+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA),
              panel.grid.minor=element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size=12),
              axis.line = element_line(colour = "black"),
              text = element_text(size = 15), #family="Source Sans Pro"),
              legend.title = element_blank(),
              legend.key = element_blank(),
        )+
        geom_text_repel(point.size=4, size=3, min.segment.length = Inf, force=0.3)+
        scale_color_manual(breaks = c("Uptake", "Release", "Unchanged"),
                           values=c("dodgerblue3", "firebrick3", "gray50"))+
        xlab("Log2fold change") + ylab("-log10(p)")+
        facet_grid(~factor(time, levels=c('pre', 'exercise', 'rec60', 'rec180')))

}


#' One sample t-tests for all time points individually
#'
#' @param df
#'
#' @return one excel results sheet for every timepoint
#' @export
#'
#' @examples
one_sample_t_tests <- function(df){

#Create a vector to loop over
timepoints <- metadata %>%
    dplyr::distinct(time) %>%
    pull(time)

for (i in timepoints){

    #Run t-tests
    results_one_sample <- {{df}} %>%
        dplyr::filter(time == i) %>%
        dplyr::group_by(time, metabolite) %>%
        t_test(av_dif_log2 ~ 1, mu=0, detailed=TRUE) %>% #One sample t-test
        dplyr::select(metabolite, time, estimate, p) %>%
        dplyr::rename(logFC=estimate,
                      P.Value = p)

    #Add extra columns
    results_one_sample <- results_one_sample %>%
        dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
        dplyr::mutate(q = qvalue::qvalue(.$P.Value, lambda=0)$qvalues) %>%
        dplyr::mutate("-log10p" = -log10(.$P.Value)) %>% #Adds -log10p column for volcano plot
        dplyr::mutate(regulated = case_when(xiao < 0.05 ~ "+")) %>%   #Adds a column called "regulated" with or without '+' depending on xiao value
        dplyr::relocate(metabolite, .before=logFC) %>%
        dplyr::arrange(desc(logFC)) %>%  #Arrange in descending order based on logFC
        dplyr::mutate(qiao = qvalue::qvalue(.$xiao, lambda=0)$qvalues)


    #Write to excel sheet
    writexl::write_xlsx(results_one_sample, here::here(paste0('data/results/results_one_sample_', i, '.xlsx')))

                                               }
}

#' One sample t-tests for all time points individually (medianscaled)
#'
#' @param df medianscaled df
#'
#' @return one excel results sheet for every timepoint (medianscaled)
#' @export
#'
#' @examples
one_sample_t_tests_medianscaled <- function(df, valid_number){

    #Create a vector to loop over
    timepoints <- metadata %>%
        dplyr::distinct(time) %>%
        pull(time)

    for (i in timepoints){

        #Make list of number of valid values
        valid_values <- {{df}} %>%
            dplyr::filter(time == i) %>%
            dplyr::group_by(peptide) %>%
            dplyr::summarise(valid_values = sum(!is.na(av_ratio))) %>%
            dplyr::filter(valid_values >= valid_number) #Filters so remaining list is of peptides with 3 or more valid values

        #Run t-tests
        results_one_sample_medianscaled <- {{df}} %>%
            dplyr::filter(time == i) %>%
            dplyr::filter(peptide %in% valid_values$peptide) %>% #filter so only peptides present in "valid_values" (from above) are retained
            dplyr::group_by(time, peptide) %>%
            t_test(av_ratio ~ 1, mu=0, detailed=TRUE) %>% #One sample t-test
            dplyr::select(peptide, time, estimate, p) %>%
            dplyr::rename(logFC=estimate,
                          P.Value = p)

        #Add extra columns
        results_one_sample_medianscaled <- results_one_sample_medianscaled %>%
            dplyr::mutate(xiao=10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
            dplyr::mutate(q = qvalue::qvalue(.$P.Value)$qvalues) %>%
            dplyr::mutate("-log10p" = -log10(.$P.Value)) %>% #Adds -log10p column for volcano plot
            dplyr::mutate(regulated = case_when(xiao < 0.05 ~ "+")) %>%   #Adds a column called "regulated" with or without '+' depending on xiao value
            dplyr::relocate(peptide, .before=logFC) %>%
            dplyr::arrange(desc(logFC)) %>%  #Arrange in descending order based on logFC
            dplyr::mutate(qiao = qvalue::qvalue(.$xiao)$qvalues)

        #Write to excel sheet
        writexl::write_xlsx(results_one_sample_medianscaled, here::here(paste0('data/results/results_one_sample_', i, '_medianscaled.xlsx')))

    }
}


#' Visualize time course of a single peptide
#'
#' @param peptide of choice. With brackets
#' @param df either df_long_ratios or df_long_ratios_medianscaled
#'
#' @return Ggplot figure
#' @export
#'
#' @examples
time_course_single <- function(df, metabolite){

    #Add numerical time column for easier plotting
    df <- {{df}} %>%
        dplyr::mutate(time_num = case_when(
            time == "pre" ~ "-70",
            time == "baseline" ~ "0",
            time == "exercise" ~ "30",
            time == "rec60" ~ "90",
            time == "rec180" ~ "210")
        ) %>%
        dplyr::mutate(time_num = as.numeric(time_num)) #change to numeric

    df %>%
        dplyr::filter(metabolite == {{metabolite}}) %>%
        ggplot(aes(x = time_num, y = av_ratio, group = 1)) +
        stat_summary(aes(x = time_num, y = av_ratio, group = id, color = "grey"),
                     fun = "mean", geom = "line", size = 1, na.rm = TRUE, show.legend = FALSE) +
        stat_summary(fun = "mean", geom = "line", size = 3, na.rm = TRUE, color = "red") +
        geom_hline(yintercept = 0, linetype = 'dashed') +
        scale_x_continuous(breaks = c(-70, 0, 30, 90, 210),
                           labels = c("Pre", "Base\nline", "Exercise", "1 h \nRecovery", "3 h \nRecovery")) +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            plot.background = element_blank(),
            axis.line = element_line(colour = "black"),
            text = element_text(size = 15),
            axis.title = element_text(size = 12),
            legend.title = element_blank(),
            legend.key = element_blank()
        ) +
        labs(y = "av_ratio", x = "Time", title = {{metabolite}})
}



#' All time points individually vs. baseline
#'
#' @param se_ratios
#'
#' @return
#' @export
#'
#' @examples
limma_analysis <- function(se_ratios, metadata) {

    timepoints <- {{metadata}} %>%
        dplyr::distinct(time) %>%
        dplyr::filter(time != "pre") %>%
        pull(time)

    all_results <- list()  # Create an empty list to store results

    for (timepoint in timepoints) {
        se_ratios_x <- se_ratios[, se_ratios$time == "baseline" | se_ratios$time == as.character(timepoint)]

        #Nedenstående subsetter ikke for de korrekte kolonner
        design <- model.matrix(~0 + se_ratios_x$timepoint)[, 1:2]
        colnames(design) <- c("baseline", as.character(timepoint))

        correlation <- duplicateCorrelation(assay(se_ratios_x), design, block = se_ratios_x$id)

        contrast <- makeContrasts(as.character(timepoint) - baseline, levels = design)

        fit <- eBayes(lmFit(assay(se_ratios), design, block = se_ratios_x$id, correlation = correlation$consensus))
        fit2 <- eBayes(contrasts.fit(fit, contrast))

        results <- topTable(fit2, coef = 1, number = Inf, sort.by = "logFC") %>%
            dplyr::mutate(xiao = 10^-(sqrt(log10(1/(P.Value^logFC))^2))) %>%
            dplyr::mutate(peptide = row.names(.),
                          q = qvalue(.$P.Value)$qvalues,
                          qiao = qvalue(.$xiao)$qvalues,
                          "-log10p" = -log10(.$P.Value),
                          regulated_xiao = ifelse(xiao < 0.05, "+", ""),
                          regulated_qiao = ifelse(qiao < 0.05, "+", ""),
                          regulated_q = ifelse(q < 0.05, "+", "")
            ) %>%
            arrange(desc(logFC))

        all_results[[as.character(timepoint)]] <- results  # Append results to the list
    }

    return(all_results)
}

