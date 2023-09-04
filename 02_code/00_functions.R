################################################################################
############################# Build object #####################################
################################################################################
## deseq_obj:
#     function that builds the dds object for the different model with
#     prefiltering, differential expression analysis
## input:
#     gse: object Large RangedSummarizedExperiment
#     model:  1 or 2; corresponding to the method of interest
#     threshold of rowsums: numerical; minimum number if counts per genes overall
#     filter: filter out genes with zero counts for more than 50% of the samples
## output:
#     DESeq obj


deseq_obj <- function(gse, model, threshold = 10, filter = T){

  if(model == 1){dds <- DESeqDataSet(gse, design = ~ diet)}
  if(model == 2){dds <- DESeqDataSet(gse, design = ~ diet + weeks + diet:weeks)}

  if(filter){
    # filter out genes with zero counts for more than 50% of the samples
    dds <- dds[which(rowSums(counts(dds) == 0)/ncol(dds) <= 0.5), ]
  } else {
    # filter out genes with less than threshold counts overall
    dds <- dds[rowSums(counts(dds)) >= threshold, ]
  }
  # Differential expression analysis
  dds <- DESeq(dds)
}
################################################################################
############################## Method I ########################################
################################################################################
# function to obtain differentially expressed gens for Method I
#
# deg: function that returns differentially expressed genes, separated by
#           up and down regulated genes
# input: dds = dds object
#        name = coefficient of interest, (resultsNames(dds))
#        pthreshold = pvalue threshold for filtering
#        log2threshold = log2foldchange threshold for filtering
#        shrinkage = T (default) or F, if shrinkage should be done or not
#
# output: results list of up and down regulated genes


de_mod1 <- function(dds, name, pthreshold = 0.05, log2threshold = log2(1.5),
                    shrinkage = T){

  if(shrinkage){
    res <- lfcShrink(dds, coef = name, type = "apeglm")
  } else {
    res <- results(dds, name = name)
  }

  res_p <- subset(res, res$padj < pthreshold)

  res_up <- subset(res_p, res_p$log2FoldChange > log2threshold)
  res_down <- subset(res_p, res_p$log2FoldChange < -log2threshold)

  return(list(up = res_up, down = res_down))
}

################################################################################
# gene_names: function to assign the gene names faster
#
# input: list of genes
#
# output: gene names of list


gene_names <- function(liste){
  gene_names <- lapply(liste, rownames)
  return(gene_names)
}

################################################################################
############################## Why filtering ###################################
################################################################################
# function for filtering
# plots a histogram for the values of the coefficient of interest
# input:
#     gse: output from deseq_obj
#     variable:  name of coefficient of interest
# ouput:
#     histogram

plot_hist <- function(gse, variable){

  df <- reshape2::melt(coef(gse))
  colnames(df) <- c("gene", "variable", "value")

  if(variable != "all"){
    df <- df[df$variable == variable,]
  }

  ggplot(data = df, aes(x= value))+
          geom_histogram(bins = 500)+
          xlim(c(-3,3))+
          theme_classic()+
          ylim(c(0,700))+
          theme(text = element_text(size = 20))+
          theme_classic()+
          theme(axis.text.x = element_text(size = 13))+
          theme(axis.text.y = element_text(size = 13))+
          theme(axis.title.x = element_text(size = 15)) +
          theme(axis.title.y = element_text(size = 15),
                plot.title = element_text(size = 15))

}
################################################################################
############################## Method II #######################################
################################################################################
# function to obtain differentially expressed gens for method ii

# input:
# x: name of the contrast of interest
# dds: DESeq object
# p_val: p-value threshold for adjusted p-values, default 0.05
# thresholdfold: for log2FoldChange, default log2(1.5)


# output:


de_mod2 <- function(x, dds, p_val = 0.05, thresholdfold = log2(1.5)){

  # shrinkage
  res_ia_x <- lfcShrink(dds, coef = x, type = "apeglm")

  # filter genes with threshold for log fold and significance
  diff_ex <- which(
    res_ia_x$padj < p_val & abs(res_ia_x$log2FoldChange) > thresholdfold)

  return(res_ia_x[diff_ex, ])

}

################################################################################
# plot for Method II all weeks

# function plots model for one gene g,
# arrows show the different parameters
# points display the observed counts,
# traingles sre mean of counts

# input:
# dds: DESeqDataSet object
# g: natural number between 1 and number of rows of DESeqData Set
# ref: which week is the reference
# point: should the observations be plottet or not
# start: lower end of intercept arrow
# s_arrow: size  arrows
# s1 : size observation points
# s2 : size mean triangle
# sl : size line

# output:
# plot


plot_it <- function(g, dds, ref, point = T, start, s_arrow  = 1.1, s1 = 3, s2 = 4, sl = 1.3){

  g1 <- "#FEC44F"
  g2 <- "#FE9929"
  g3 <- "#EC7014"
  g4 <- "#CC4C02"
  g5 <- "#993404"
  b1 <- "#9ECAE1"
  b2 <- "#6BAED6"
  b3 <- "#4292C6"
  b4 <- "#2171B5"
  b5 <- "#08519C"
  b0 <- "yellowgreen"
  a  <- "gray25"
  l1 <- "gray56"
  l2 <- "antiquewhite3"


  # data.frame of log2 counts
  log_count <- data.frame(y = log2(counts(dds, normalize = T)[g,]),
                          x = as.numeric(sub("w.*", "",
                                             names(counts(dds)[g,]))),
                          diet =  as.factor(gsub(".*w-(.+)-M.*", "\\1",
                                                 names(counts(dds)[g,]))))
  # change level HFD to WD
  levels(log_count$diet)[levels(log_count$diet) == "HFD"] = "WD"

  # create column for plotting containing week and diet together
  log_count$index <- paste0(log_count$x, log_count$diet)

  # mean log2 count for each time point per diet
  mean_data <- data.frame(x = rep(unique(log_count$x), each = 2),
                          meanCount = sapply(unique(log_count$index),
                                             function(i)
                                               log2(mean(2^log_count[log_count$index == i,]$y))),
                          diet = as.factor(rep(c("WD", "SD"), 6)))




   # combine all data in one data.frame
  all_data <- data.frame(x = c(log_count$x, mean_data$x),
                         y = c(log_count$y, mean_data$meanCount),
                         diet = c(log_count$diet, mean_data$diet),
                         count = as.factor(
                           c(rep("observed", nrow(log_count)),
                             rep("mean", nrow(mean_data)))))

  all_data$count <- factor(all_data$count, levels = c("observed", "mean"))

  # data.frame containing the values for start and end points of the arrows
  df <- data.frame(
    # weeks with reference week in first position (x-axis values):
    x = rep(c(ref,
              c(3, 6, 30, 36, 42, 48)[- which(c(3, 6, 30, 36, 42, 48) == ref)]),
            3),


    # response (y-axis):
    y = c(
      # diet sd for all weeks (for all but the reference week add
      # intercept):
      coef(dds)[g,1],
      sapply(3:7, function(x) coef(dds)[g, x] + coef(dds)[g,1]),

      # intercept and standarddiet + effect of WD
      coef(dds)[g,1] + coef(dds)[g,2],
      sapply(3:7, function(x) coef(dds)[g,x] + coef(dds)[g,2] + coef(dds)[g,1]),

      # WD values with interaction
      coef(dds)[g,1] + coef(dds)[g,2],
      sapply(3:7, function(x) coef(dds)[g, x] + coef(dds)[g, 2] +
               coef(dds)[g, 1] + coef(dds)[g, x + 5])
    ),

    # level for colors
    diet = as.factor(c(rep("SD", 6), rep("sd_wd", 6),
                       rep("WD", 6)))
  )

  # for the line for the diets only WD and sd is needed
  df_line <-  df[df$diet != "sd_wd", ]
  df_line$diet <-  droplevels(df_line$diet)


  #data.frame for the week effect arrows
  df_week <- data.frame(x = c(3, 6, 30, 36, 42, 48)[- which(c(3, 6, 30, 36, 42,
                                                              48) == ref)],
                        xend = c(3, 6, 30, 36, 42,
                                 48)[- which(c(3, 6, 30, 36, 42, 48) == ref)],
                        yend = df[2:6, 2],
                        y = rep(coef(dds)[g, 1], 5))

  # data for the interaction arrows
  df_int_arrow <- data.frame(x = c(3, 6, 30, 36, 42,
                                   48)[- which(c(3, 6, 30, 36, 42, 48) == ref)],
                             xend = c(3, 6, 30, 36, 42,
                                      48)[- which(c(3, 6, 30, 36, 42, 48) == ref)],
                             y = df[8:12, 2],
                             yend = df[14:18,2])

  # main effect arrows
  df_hfd_ef <- data.frame(x = c(ref, c(3, 6, 30, 36, 42,
                                       48)[- which(c(3, 6, 30, 36, 42, 48) == ref)]),
                          xend = c(ref, c(3, 6, 30, 36, 42,
                                          48)[- which(c(3, 6, 30, 36, 42, 48) == ref)]),
                          yend = df[7:12, 2],
                          y = df[1:6,2])

  ggplot(df, aes(x = x, y = y)) +
    theme_classic() +
    scale_color_manual(values = c(l1, l2))+
    # lines for sd and hfd
    geom_line(data = df_line, aes(x = x, y = y, col = diet), size = sl) +
    # plot observations:
    {if(point == T){
      geom_jitter(data = all_data, aes(x = x, y = y, col = diet, shape = count),
                  alpha = c(rep(1, nrow(log_count)), rep(1, nrow(mean_data))),
                  size = c(rep(s1, nrow(log_count)), rep(s2, nrow(mean_data))),
                  width = 0.5)
    }}+

    # arrow for intercept
    {if(coef(dds)[g,1] <= 0){
      geom_segment(aes(x = ref + 0.6, y = 0, xend = ref + 0.6,
                       yend = coef(dds)[g,1]),
                   arrow = arrow(length = unit(0.2, "cm")), col = b0,
                   size = s_arrow)}
      else{
        geom_segment(aes(x = ref + 0.6, y = start,
                         xend = ref + 0.6,
                         yend = coef(dds)[g,1]),
                     arrow = arrow(length = unit(0.2, "cm")), col = b0,
                     size = s_arrow)}}+

    # arrow for week effect
    geom_segment(data = df_week, aes(x = x - 0.7, xend = xend - 0.7, y = y,
                                     yend = yend),
                 arrow = arrow(length = unit(0.2, "cm")),
                 col = c(b1, b2, b3  ,b4, b5), size = s_arrow) +


    # arrow hfd effect
    geom_segment(data = df_hfd_ef, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.2, "cm")), col = a,
                 size = s_arrow) +


    # hfd with interaction
    geom_segment(data = df_int_arrow, aes(x = x + 0.7, xend = xend + 0.7, y = y,
                                          yend = yend),
                 arrow = arrow(length = unit(0.2, "cm")),
                 col = c(g1, g2, g3 , g4 , g5), size = s_arrow) +


    labs(x = "weeks", y = "log2(normalized count)") +
    scale_x_continuous(limits = c(0, 53), breaks = seq(0, 50, by = 5))+
    theme(axis.text.x = element_text(size = 15))+
    theme(axis.text.y = element_text(size = 15))+
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16))+
    theme(legend.position="bottom")+
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_text(colour="black", size = 15,
                                      face="bold"))+
    theme(legend.text = element_text(colour="black", size = 15))
}

################################################################################
# function to plot different forms of interaction

# input:
#   df: data frame contatining one column for x values, one for y values, and one for color
#   plot_title: title for the plot
# output: plot

scheme_interaction_plot <- function(df, plot_title){
  ggplot(df, aes(x = x, y = y, col = col)) +
    geom_point(col = c( "blue",  "blue", "red", "red"), size= 2) +
    geom_line(aes( x = x, y = y, colour = col),  size = 1.5)+
    theme_classic()+
    scale_color_manual("", breaks = c("Group 0","Group 1"), values = c( "blue","red"))+
    theme(axis.text.x = element_text(size = 22))+
    theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())+
    theme(axis.title.x = element_text(size = 24)) +
    theme(legend.key.size = unit(3, "cm"), legend.text = element_text(size=20))+
    theme(axis.title.y = element_text(size = 24))+
    ggtitle(plot_title) +
    ylab("")+
    scale_x_discrete(limits = c("low", "high"), name = "")+
    theme(plot.title = element_text(size = 25))
}
################################################################################
# function to plot different possibilities for significance between Method I and II
# generates normally distributed data

# input:
# y1, y2, y3, y4: different y-values to create scenarios with different interaction effects
# v: variance for generated data
# caption1, caption2: caption for plots, used for either "DEG" or "not DEG"
# l1: linetype for model 1 week 6 diet effect
# l2: linetype for model 1 week 3 diet effect
# l2: linetype for model 3 interaction effect
# l4: linetype for model 3 main effect
# l1-l3: "dashed" or "solid"
# number: annotate plot to corresponding area in scatter plot (1-7)

scheme_model <- function(y1, y2, y3, y4, v = 0.1, caption1, caption2, l1, l2, l3,l4,  number){
  # values for the plot
  arrow_size = 1.5
  point_size = 2
  # x values for the arrows:
  x1 = 0.92
  x2 = 1.42
  x3 = 0.95
  x4 = 1.45
  x5 = 1.39

  # create data
  set.seed(3)
  data_scheme <- data.frame(
    "y" = c(mvrnorm(5, y1, v, empirical = T), mvrnorm(5, y2, v, empirical = T),
            mvrnorm(5, y3, v, empirical = T), mvrnorm(5, y4, v, empirical = T)),
    "week" = (c(rep(1, 5), rep(1.07, 5), rep(1.5, 5), rep(1.57, 5))),
    "group" = as.factor(c(rep(c(rep("SD", 5), rep("WD", 5)), 2))))

  # plot mehtod i
  pl1 <- ggplot(data_scheme, aes(x = week, y = y, col = group)) +
    geom_point(size = point_size) +
    theme_classic() +
    #intercept week3
    geom_segment(aes(x = x1, xend = x1, y = 0, yend = y1),
                 arrow = arrow(length = unit(0.2, "cm")), col = "black", size = arrow_size)+
    # intercept week6
    geom_segment(aes(x = x2, xend = x2, y = 0, yend = y3),
                 arrow = arrow(length = unit(0.2, "cm")), col = "black", size = arrow_size)+
    # main effect week 3
    geom_segment(aes(x = x3, xend = x3, y = y1, yend = y2),
                 arrow = arrow(length = unit(0.2, "cm")), col = "darkgreen",
                 size = arrow_size, linetype = 1)+
    {if(l1 == F){
      geom_segment(aes(x = x3, xend = x3, y = y1, yend = y2 - 0.03),
                   col = "white",
                   size = arrow_size, linetype = 9)
    }} +
    # main effect week6
    geom_segment(aes(x = x4, xend = x4, y = y3, yend = y4),
                 arrow = arrow(length = unit(0.2, "cm")), col = "darkgreen",
                 size = arrow_size, linetype = 1)+
    {if(l2 == F){
      geom_segment(aes(x = x4, xend = x4, y = y3  , yend = y4 + 0.04),
                   col = "white",
                   size = arrow_size, linetype = 9)
    }} +

    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),axis.ticks.y=element_blank()) +
    # caption weather DEG or not
    {if(caption1 != "not DEG"){
      theme(plot.caption = element_text(color = "seagreen", size = 12, face = "bold"))}
      else{
        theme(plot.caption = element_text(color = "red3", size = 12, face = "bold"))}}+
    labs(caption = caption1) +
    theme(plot.caption=element_text(size=20))+
    theme(legend.title=element_blank(), legend.text=element_text(size=20)) +
    ggtitle(number)+
    theme(plot.title = element_text(size = 20)) +
    theme(legend.position="none") +
    theme(plot.background=element_rect(color="black"))

  # plot for method ii
  pl10 <- ggplot(data_scheme, aes(x = week, y = y, col = group)) +
    geom_point(size = point_size) +
    theme_classic()  +
    # intercept
    geom_segment(aes(x = x1, xend = x1, y = 0, yend = y1),
                 arrow = arrow(length = unit(0.2, "cm")), col = "black", size = arrow_size)+
    # week effect
    geom_segment(aes(x = x5, xend = x5, y = y1, yend = y3),
                 arrow = arrow(length = unit(0.2, "cm")), col = "azure4", size = arrow_size)+
    # main effect week 3
    geom_segment(aes(x = x3, xend = x3, y = y1, yend = y2),
                 arrow = arrow(length = unit(0.2, "cm")), col = "dodgerblue3",
                 size = arrow_size, linetype = 1)+
    {if(l3 == F){
      geom_segment(aes(x = x3, xend = x3, y = y1, yend = y2 - 0.03),
                   col = "white", size = arrow_size, linetype = 9)
    }}+
    # main effect week 6
    geom_segment(aes(x = x2, xend = x2, y = y3, yend = y3 + y2 - y1 ),
                 arrow = arrow(length = unit(0.2, "cm")), col = "dodgerblue3",
                 size = arrow_size, linetype = 1) +
    {if(l3 == F){
      geom_segment(aes(x = x2, xend = x2, y = y3, yend = y3 + y2 - y1 - 0.03),
                   col = "white",
                   size = arrow_size, linetype = 9)
    }} +
    # interaction effect
    geom_segment(aes(x = x4, xend = x4, y = y3 + y2 - y1, yend = y4),
                 arrow = arrow(length = unit(0.2, "cm")), col = "orchid3",
                 size = arrow_size, linetype = 1) +
    {if(l4 == F){
      geom_segment(aes(x = x4, xend = x4, y = y3 + y2 - y1, yend = y4  - 0.06),
                   col = "white",
                   size = arrow_size, linetype = 9)
    }} +
    # line for week 3 reference
    geom_segment(aes(x = x1, xend = 1.6, y = y1, yend = y1), col = "grey",
                 linetype = "dashed")+
    theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),
          axis.title.x=element_blank()) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y=element_blank()) +

    {if(caption2 != "not DEG"){
      theme(plot.caption = element_text(color = "seagreen", size = 12, face = "bold"))}
      else{
        theme(plot.caption = element_text(color = "red3", size = 12, face = "bold"))} }+

    labs(caption = caption2) +
    theme(plot.caption=element_text(size=20))+
    theme(legend.title=element_blank(), legend.text=element_text(size=20)) +
    ggtitle(" ") +
    theme(plot.title = element_text(size = 20))+
    theme(legend.position="none")+
    theme(plot.background=element_rect(color="black"))

  return(list(pl1, pl10))
}

################################################################################
# plot for Method II only week 3 and 6

# function plots model for one gene g,
# arrows show the different parameters
# points display the observed counts
# triangle, mean counts

# input:
# dds: DESeqDataSet object
# g: natural number between 1 and number of rows of DESeqData Set
# ref: which week is the reference
# point: should the observations be plottet or not
# start: lower end of intercept arrow
# s_arrow: size  arrows
# s1 : size observation points
# s2 : size mean triangle
# sl : size line
# lim1: lower limit for y axis
# lim2: upper limit for y axis
# main: true false, if main effect is bigger log2(1.5) -> gestrichelt
# interac: true false, if interac effect is bigger log2(1.5) -> gestrichelt

# output:
# plot


plot_it_3_6 <- function(g, dds, ref, point = T, start, s_arrow  = 2.5, s1 = 4, s2 = 5, sl = 1.3, lim1, lim2, main = T, interac = T){


    # data.frame of log2 counts
    log_count <- data.frame(y = log2(counts(dds, normalize = T)[g,]),
                            x = as.numeric(sub("w.*", "",
                                               names(counts(dds)[g,]))),
                            diet =  as.factor(gsub(".*w-(.+)-M.*", "\\1",
                                                   names(counts(dds)[g,]))))[1:22,]

    # change level HFD to WD
    levels(log_count$diet)[levels(log_count$diet) == "HFD"] = "WD"

    # create column for plotting containing week and diet together
    log_count$index <- paste0(log_count$x, log_count$diet)[1:22]

    # mean log2 count for each time point per diet
    mean_data <- data.frame(x = rep(unique(log_count$x), each = 2),
                            meanCount = sapply(unique(log_count$index),
                                               function(i)
                                                 log2(mean(2^log_count[log_count$index == i,]$y))),
                            diet = as.factor(rep(c("WD", "SD"), 2)))

    # combine all data in one data.frame
    all_data <- data.frame(x = c(log_count$x, mean_data$x),
                           y = c(log_count$y, mean_data$meanCount),
                           diet = c(log_count$diet, mean_data$diet),
                           count = as.factor(
                             c(rep("observed", nrow(log_count)),
                               rep("mean", nrow(mean_data)))))

    all_data$count <- factor(all_data$count, levels = c("observed", "mean"))

    # data.frame containing the values for start and end points of the arrows


    #data.frame for the week effect arrows
    df_week <- data.frame(x = 6,
                          xend = 6,
                          yend = coef(dds)[g, 1] + coef(dds)[g, 3],
                          y = coef(dds)[g, 1])

    # main effect arrows
    df_hfd_ef <- data.frame(x = c(3, 6),
                            xend = c(3, 6),
                            yend = c(coef(dds)[g, 1] + coef(dds)[g, 2],
                                     coef(dds)[g, 1] + coef(dds)[g, 3] + coef(dds)[g, 2]),
                            y = c(coef(dds)[g, 1], coef(dds)[g, 1] + coef(dds)[g, 3]))


    # data for the interaction arrows
    df_int_arrow <- data.frame(x = 6,
                               xend = 6,
                               y = coef(dds)[g, 1] + coef(dds)[g, 3] + coef(dds)[g, 2],
                               yend = coef(dds)[g, 1] + coef(dds)[g, 3] + coef(dds)[g, 2] + coef(dds)[g, 8])


    # values to shift arrows -> arrows next to each other not on top
    x1 = 1.5
    x2 = 1
    x3 = 0.5

    # start plot
    ggplot(all_data, aes(x = x, y = y)) +
      theme_classic() +
      # plot observations:
      geom_jitter(data = all_data[all_data$count != "mean",],
                  aes(x = x, y = y, col = diet, shape = count),
                  alpha = c(rep(0.75, nrow(log_count))),
                  shape = c(rep(16, nrow(log_count))),
                  size = c(rep(s1, nrow(log_count))),
                  width = 0.2, height = 0) +

      geom_jitter(data = all_data[all_data$count == "mean",],
                  aes(x = x, y = y, col = diet, shape = count),
                  alpha = rep(0.9, nrow(mean_data)),
                  shape = rep(15, nrow(mean_data)),
                  size = rep(s2, nrow(mean_data)),
                  width = 0.001, height = 0) +


      # arrow for intercept, if: intercept <= 0, else: intercept > 0
      {if(coef(dds)[g,1] <= 0){
        geom_segment(aes(x = 3 - x1, y = 0, xend = 3 - x1,
                         yend = coef(dds)[g,1]),
                     arrow = arrow(length = unit(0.2, "cm")), col = "black",
                     size = s_arrow)}
        else{
          geom_segment(aes(x = 3 - x1, y = start,
                           xend = ref - x1,
                           yend = coef(dds)[g,1]),
                       arrow = arrow(length = unit(0.2, "cm")), col = "black",
                       size = s_arrow)}}+

      # arrow for week effect
      geom_segment(data = df_week, aes(x = x - x1, xend = xend - x1, y = y,
                                       yend = yend),
                   arrow = arrow(length = unit(0.2, "cm")),
                   col = "azure4", size = s_arrow) +

      # arrow hfd effect
      geom_segment(data = df_hfd_ef, aes(x = x - x2, xend = xend - x2, y = y, yend = yend),
                   arrow = arrow(length = unit(0.2, "cm")), col = "dodgerblue3",
                   size = s_arrow) +

      # if  effect < log2(1.5) line is dotted
      {if(main != T){
        geom_segment(data = df_hfd_ef, aes(x = x - x2, xend = xend - x2, y = y, yend = yend - main) ,
                     col = "white", size = s_arrow, linetype = 9, alpha = 0.8)
      }}+

      # hfd with interaction
      geom_segment(data = df_int_arrow, aes(x = x - x3, xend = xend - x3, y = y,
                                            yend = yend),
                   arrow = arrow(length = unit(0.2, "cm")),
                   col = "orchid3", size = s_arrow) +

      # if  effect < log2(1.5) line is dotted
      {if(interac != T){
        geom_segment(data = df_int_arrow, aes(x = x - x3, xend = xend - x3, y = y,
                                              yend = yend - interac),
                     col = "white", size = s_arrow, linetype = 9, alpha = 0.8)
      }}+

      # dashed horizontal line for week effect starting point
      geom_segment(aes(x = 3 -x1, xend = 6 -x1, y = coef(dds)[g,1], yend = coef(dds)[g,1]), col = "grey",  linetype = "dashed")+

      labs(x = "", y = "")+
      scale_x_continuous(limits = c(1.5, 6.5), breaks = c(3, 6))+
      scale_y_continuous(limits = c(lim1, lim2))+
      ggtitle("")+
      theme(axis.text.x = element_text(size = 15))+
      theme(axis.text.y = element_text(size = 15))+
      theme(axis.title.x = element_text(size = 16)) +
      theme(axis.title.y = element_text(size = 16))+
      theme(legend.position="bottom")+
      guides(shape = guide_legend(override.aes = list(size = 3))) +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      theme(legend.title = element_text(colour="black", size = 15,
                                        face="bold"))+
      theme(plot.background=element_rect(color="black"))+
      theme(legend.position="none")+
      theme(legend.text = element_text(colour="black", size = 15))
}

################################################################################
# plot for model 3 only week 3 and 6

# function plots model for one gene g,
# arrows show the different parameters
# points display the observed counts

# input:
# dds: DESeqDataSet object from method ii
# week3: DESeqDataSet object for week 3 model 1
# week6: DESeqDataSet object for week 6 model 1
# g: natural number between 1 and number of rows of DESeqData Set
# ref: which week is the reference
# point: should the observations be plottet or not
# start: lower end of intercept arrow
# s_arrow: size  arrows
# s1 : size observation points
# s2 : size mean triangle
# sl : size line
# lim1: lower limit for y axis
# lim2: upper limit for y axis
# main3: true false, if main effect in week3 is bigger log2(1.5) -> dotted
# main6: true false, if main effect in week 6 is bigger log2(1.5) -> dotted



# output:
# plot


plot_it_3_6_mod1 <- function(g, week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]],
                             dds, ref, point = T, start, s_arrow  = 2.5, s1 = 4, s2 = 5,
                             sl = 1.3, lim1, lim2, main3 = T, main6 = T, number){

  # data.frame of log2 counts
  log_count <- data.frame(y = c(log2(counts(week3, normalize = T)[g,]),
                                log2(counts(week6, normalize = T)[g,])),
                          x = c(as.numeric(sub("w.*", "",
                                               names(counts(week3)[g,]))),
                                as.numeric(sub("w.*", "",
                                               names(counts(week6)[g,])))),
                          diet =  c(as.factor(gsub(".*w-(.+)-M.*", "\\1",
                                                   names(counts(week3)[g,]))),
                                    as.factor(gsub(".*w-(.+)-M.*", "\\1",
                                                   names(counts(week6)[g,])))))


  # change level HFD to WD
  levels(log_count$diet)[levels(log_count$diet) == "HFD"] = "WD"
  # create column for plotting containing week and diet together
  log_count$index <- paste0(log_count$x, log_count$diet)

  # mean log2 count for each time point per diet
  mean_data <- data.frame(x = rep(unique(log_count$x), each = 2),
                          meanCount = sapply(unique(log_count$index),
                                             function(i)
                                               log2( mean(2^log_count[log_count$index == i,]$y))),
                          diet = as.factor(rep(c("WD", "SD"), 2)))

  # combine all data in one data.frame
  all_data <- data.frame(x = c(log_count$x, mean_data$x),
                         y = c(log_count$y, mean_data$meanCount),
                         diet = c(log_count$diet, mean_data$diet),
                         count = as.factor(
                           c(rep("observed", nrow(log_count)),
                             rep("mean", nrow(mean_data)))))

  all_data$count <- factor(all_data$count, levels = c("observed", "mean"))


  # main effect arrows
  df_hfd_ef <- data.frame(x = c(3, 6),
                          xend = c(3, 6),
                          yend = c(coef(week3)[g, 1] + coef(week3)[g, 2],
                                   coef(week6)[g, 1] + coef(week6)[g, 2]),
                          y = c(coef(week3)[g, 1], coef(week6)[g, 1]))

  # values to shift arrows -> arrows next to each other not on top
  x2 = 1
  x3 = 0.5


  ggplot(all_data, aes(x = x, y = y)) +
    theme_classic() +
    # plot observations:
    geom_jitter(data = all_data[all_data$count != "mean",],
                aes(x = x, y = y, col = diet, shape = count),
                alpha = c(rep(0.75, nrow(log_count))),
                shape = c(rep(16, nrow(log_count))),
                size = c(rep(s1, nrow(log_count))),
                width = 0.2, height = 0) +

    geom_jitter(data = all_data[all_data$count == "mean",],
                aes(x = x, y = y, col = diet, shape = count),
                alpha = rep(0.9, nrow(mean_data)),
                shape = rep(15, nrow(mean_data)),
                size = rep(s2, nrow(mean_data)),
                width = 0.001, height =0) +


    # arrow for intercept week 3, if: intercept <= 0, else > 0
    {if(coef(week3)[g,1] <= 0){
      geom_segment(aes(x = 3 - x2, y = 0, xend = 3 - x2,
                       yend = coef(week3)[g,1]),
                   arrow = arrow(length = unit(0.2, "cm")), col = "black",
                   size = s_arrow) }
      else{
        geom_segment(aes(x = 3 - x2, y = start,
                         xend = 3 - x2,
                         yend = coef(week3)[g,1]),
                     arrow = arrow(length = unit(0.2, "cm")), alpha = 0.8,  col = "black",
                     size = s_arrow)}}+
    # arrow for intercept week 6, if: intercept <= 0, else > 0
    {if(coef(week6)[g,1] <= 0){
      geom_segment(aes(x = 6 - x2, y = 0, xend = 6 - x2,
                       yend = coef(week6)[g,1]),
                   arrow = arrow(length = unit(0.2, "cm")), col = "black",
                   size = s_arrow) }
      else{
        geom_segment(aes(x = 6 - x2, y = start,
                         xend = 6 - x2,
                         yend = coef(week6)[g,1]),
                     arrow = arrow(length = unit(0.2, "cm")), col = "black",
                     size = s_arrow, alpha = 0.8)}}+


    # arrow hfd effect
    geom_segment(data = df_hfd_ef, aes(x = x - x3, xend = xend - x3, y = y, yend = yend),
                 arrow = arrow(length = unit(0.2, "cm")), col =  "darkgreen",
                 size = s_arrow) +

    # when effect abs < log2(1.5) -> dotted
    {if(main3 != T){
      geom_segment(data = df_hfd_ef[df_hfd_ef$x == 3, ],
                   aes(x = x - x3, xend = xend - x3, y = y,
                       yend = yend - main3),
                   col = "white", size = s_arrow, linetype = 9, alpha = 0.8)
    }}+

    # when effect abs < log2(1.5) -> dotted
    {if(main6 != T){
      geom_segment(data = df_hfd_ef[df_hfd_ef$x == 6, ],
                   aes(x = x - x3, xend = xend - x3, y = y, yend = yend - main6),
                   col = "white", size = s_arrow, linetype = 9, alpha = 0.8)
    }}+


    labs(x = "", y = "")+
    scale_x_continuous(limits = c(2, 7), breaks = c(3, 6))+
    scale_y_continuous(limits = c(lim1, lim2))+
    ggtitle(number)+
    theme(axis.text.x = element_text(size = 15))+
    theme(axis.text.y = element_text(size = 15),
          plot.title =   element_text(size = 18))+
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16))+
    theme(legend.position="bottom")+
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_text(colour="black", size = 15,
                                      face="bold"))+
    theme(plot.background=element_rect(color="black"))+
    theme(legend.position="none")+
    theme(legend.text = element_text(colour="black", size = 15))
}

################################################################################
# plot for method i

# function plots model for one gene g,
# arrows show the different parameters
# points display the observed counts

# input:
# dds: DESeqDataSet object for model 1
# g: natural number between 1 and number of rows of DESeqData Set
# ref: which week is the reference
# point: should the observations be plottet or not
# start: lower end of intercept arrow
# s_arrow: size  arrows
# s1 : size observation points
# s2 : size mean triangle
# sl : size line

# output:
# plot


plot_no_int <- function(g, dds, ref, point = T, start, s_arrow  = 1.1, s1 = 3, s2 = 4, sl = 1.3){


  b1 <- "#9ECAE1"
  b2 <- "#6BAED6"
  b3 <- "#4292C6"
  b4 <- "#2171B5"
  b5 <- "#08519C"
  b0 <- "yellowgreen"
  a <- "gray25"
  l1 <- "gray56"
  l2 <- "antiquewhite3"

  # data.frame of log2 counts
  log_count <- data.frame(y = log2(counts(dds, normalize = T)[g,]),
                          x = as.numeric(sub("w.*", "",
                                             names(counts(dds)[g,]))),
                          diet =  as.factor(gsub(".*w-(.+)-M.*", "\\1",
                                                 names(counts(dds)[g,]))))

  # change level HFD to WD
  levels(log_count$diet)[levels(log_count$diet) == "HFD"] = "WD"
  # create column for plotting containing week and diet together
  log_count$index <- paste0(log_count$x, log_count$diet)

  # mean log2 count for each time point per diet
  mean_data <- data.frame(x = rep(unique(log_count$x), each = 2),
                          meanCount = sapply(unique(log_count$index),
                                             function(i)
                                               log2(mean(2^log_count[log_count$index == i,]$y))),
                          diet = as.factor(rep(c("WD", "SD"), 6)))



  # combine all data in one data.frame
  all_data <- data.frame(x = c(log_count$x, mean_data$x),
                         y = c(log_count$y, mean_data$meanCount),
                         diet = c(log_count$diet, mean_data$diet),
                         count = as.factor(
                           c(rep("observed", nrow(log_count)),
                             rep("mean", nrow(mean_data)))))

  all_data$count <- factor(all_data$count, levels = c("observed", "mean"))

  # data.frame containing the values for start and end points of the arrows
  df <- data.frame(
    # weeks with reference week in first position (x-axis values):
    x = rep(c(ref,
              c(3, 6, 30, 36, 42, 48)[- which(c(3, 6, 30, 36, 42, 48) == ref)]),
            2),


    # response (y-axis):
    y = c(
      # diet sd for all weeks (for all but the reference week add
      # intercept):
      coef(dds)[g,1],
      sapply(3:7, function(x) coef(dds)[g, x] + coef(dds)[g,1]),

      # intercept and standarddiet + effect of WD
      coef(dds)[g,1] + coef(dds)[g,2],
      sapply(3:7, function(x) coef(dds)[g,x] + coef(dds)[g,2] + coef(dds)[g,1])
    ),

    # level for colors
    diet = as.factor(c(rep("SD", 6),
                       rep("WD", 6)))
  )



  #data.frame for the week effect arrows
  df_week <- data.frame(x = c(3, 6, 30, 36, 42, 48)[- which(c(3, 6, 30, 36, 42,
                                                              48) == ref)],
                        xend = c(3, 6, 30, 36, 42,
                                 48)[- which(c(3, 6, 30, 36, 42, 48) == ref)],
                        yend = df[2:6, 2],
                        y = rep(coef(dds)[g, 1], 5))


  # main effect arrows
  df_hfd_ef <- data.frame(x = c(ref, c(3, 6, 30, 36, 42,
                                       48)[- which(c(3, 6, 30, 36, 42, 48) == ref)]),
                          xend = c(ref, c(3, 6, 30, 36, 42,
                                          48)[- which(c(3, 6, 30, 36, 42, 48) == ref)]),
                          yend = df[7:12, 2],
                          y = df[1:6,2])



  ggplot(df, aes(x = x, y = y)) +
    theme_classic() +
    scale_color_manual(values = c(l1, l2))+
    # lines for sd and WD
    geom_line(data = df, aes(x = x, y = y, col = diet), size = sl) +
    # plot observations:
    {if(point == T){
      geom_jitter(data = all_data, aes(x = x, y = y, col = diet, shape = count),
                  alpha = c(rep(1, nrow(log_count)), rep(1, nrow(mean_data))),
                  size = c(rep(s1, nrow(log_count)), rep(s2, nrow(mean_data))),
                  width = 0.4)
    }}+

    # arrow for intercept
    {if(coef(dds)[g,1] <= 0){
      geom_segment(aes(x = ref + 0.6, y = 0, xend = ref + 0.6,
                       yend = coef(dds)[g,1]),
                   arrow = arrow(length = unit(0.2, "cm")), col = b0,
                   size = s_arrow)}
      else{
        geom_segment(aes(x = ref + 0.6, y = floor(min(log_count$y)- 0.1),
                         xend = ref + 0.6,
                         yend = coef(dds)[g,1]),
                     arrow = arrow(length = unit(0.2, "cm")), col = b0,
                     size = s_arrow)}}+

    # arrow for week effect
    geom_segment(data = df_week, aes(x = x - 0.7, xend = xend - 0.7, y = y,
                                     yend = yend),
                 arrow = arrow(length = unit(0.2, "cm")),
                 col = c(b1, b2, b3  ,b4, b5), size = s_arrow) +


    # arrow WD effect
    geom_segment(data = df_hfd_ef, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = arrow(length = unit(0.2, "cm")), col = a,
                 size = s_arrow) +



    labs(x = "weeks", y = "log2(normalized count)") +
    scale_x_continuous(limits = c(0, 53), breaks = seq(0, 50, by = 5))+

    theme(axis.text.x = element_text(size = 15))+
    theme(axis.text.y = element_text(size = 15))+
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16))+
    theme(legend.position="bottom")+
    guides(shape = guide_legend(override.aes = list(size = 3))) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_text(colour="black", size = 15,
                                      face="bold"))+
    theme(legend.text = element_text(colour="black", size = 15))
}


################################################################################

# Functions for GO-analysis (unsing package topGO)

################################################################################

# Helper function
topGOdataObj <- function(interesting.genes, gene.names = gene_id_all, node.size = 5) {

  gene_list <- factor(as.integer(gene.names %in% interesting.genes))
  names(gene_list) <- gene.names

  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = gene_list,
                nodeSize = node.size,
                annot = annFUN.org,
                mapping ="org.Mm.eg.db",
                ID = "ensembl")
}


################################################################################
# function for go analysis for groups from comparison plot
# input:
#   intersting: name of interesting genes
#   all_genes: name of all genes
#   methods: method used for the runTest funtion
#   topNodes: natural number, how many Nodes are of interest


# output:
#   data frame containing  the go results with adjusted p-values

topGO_analysis <- function(interesting, all_genes, methods = "elim", topNodes = 1000){



  # create topGo object
  go_data <- topGOdataObj(interesting.genes = interesting,
                          gene.names = all_genes)

  # apply method to topGo Object
  res <- sapply(methods,
                function(x) runTest(go_data, algorithm = x, statistic = "fisher"))

  # res$weight.log <- runTest(GOdata, algorithm = "weight", statistic = "fisher", sigRatio = "log")
  # res$parentchild.intersect <- runTest(GOdata, algorithm = "parentchild", statistic = "fisher", joinFun = "intersect")
  
  
  if(methods == "elim"){
    res$gen_table <- GenTable(go_data,
                              Fisher.elim = res$elim,
                              orderBy = "Fisher_elim", topNodes = topNodes)
    
    
    res$gen_table$Fisher.elim <- as.numeric(res$gen_table$Fisher.elim)
    if(any(is.na(res$gen_table$Fisher.elim))){
      res$gen_table$Fisher.elim[is.na(res$gen_table$Fisher.elim)] <- 1e-30
    }
    ## adjust the p-values
    res$gen_table$Fisher.elim_adjust <- p.adjust(res$gen_table$Fisher.elim, method = "fdr")
  }
  
  
  if(methods == "classic"){
    res$gen_table <- GenTable(go_data,
                              Fisher.elim = res$classic,
                              orderBy = "Fisher_elim", topNodes = topNodes)
    
    
    res$gen_table$Fisher.elim <- as.numeric(res$gen_table$Fisher.elim)
    if(any(is.na(res$gen_table$Fisher.elim))){
      res$gen_table$Fisher.elim[is.na(res$gen_table$Fisher.elim)] <- 1e-30
    }
    ## adjust the p-values
    res$gen_table$Fisher.elim_adjust <- p.adjust(res$gen_table$Fisher.elim, method = "fdr")
  }
  return(res)
  
}

###############################################################################

# Dot plot for Results from go analysis
# input:
#     - dfgene: value 'gen_table' from output of topGO_analysis2
#     - n: natural number of top n groups with smallest adjusted p-value
#     - title: optional; title for the plot

# output: dotplot (ggplot-object)


dot_plot <- function(dfgene, n, add_title=""){
  df <- dfgene[1:n,]

  if(sum(grepl("antigen processing and presentation of e...", df$Term)) > 1){
    df$Term[grep("antigen processing and presentation of e...",
                 df$Term)] <-
      c("antigen processing/presentation via MHC class I", #antigen processing and presentation of exogenous peptide antigen via MHC class I
        "antigen processing/presentation via MHC class II") #antigen processing and presentation of exogenous peptide antigen via MHC class II
  }
  df$Term <- gsub("positive", "pos.", df$Term)
  df$Term <- gsub("negative", "neg.", df$Term)
  df$Term <- gsub("regulation", "reg.", df$Term)
  df$Term <- gsub("necrosis fa...", "necrosis factor", df$Term)
  df$Term <- gsub("cas...", "cascade", df$Term)
  df$Term <- gsub("intrinsic apoptotic signal...", "intrinsic apoptotic signaling pathway", df$Term)
  df$Term <- gsub("interferon-gamma ...", "interferon-gamma production", df$Term)
  df$Term <- gsub("interleukin-6 pro...", "interleukin-6", df$Term)
  df$Term <- gsub("interleukin-1 bet...", "interleukin-1", df$Term)
  df$Term <- gsub("defense response", "defense resp.", df$Term)
  df$Term <- gsub("bacter\\.\\.\\.", "bacterium", df$Term)
  df$Term <- gsub("inflammatory resp\\.\\.\\.", "inflammatory response", df$Term)
  df$Term <- gsub("pat\\.\\.\\.|proc\\.\\.\\.", "", df$Term)
  df$Term <- gsub("transport|transp...", "transp.", df$Term)
  df$Term <- gsub("transmis\\.\\.\\.", "transmission", df$Term)
  df$Term <- gsub("e1\\.\\.\\.", "e1", df$Term)
  df$Term <- gsub("e2\\.\\.\\.", "e2", df$Term)
  df$Term <- gsub("act\\.\\.\\.", "activity", df$Term)
  df$Term <- gsub("he\\.\\.\\.", "healing", df$Term)
  df$Term <- gsub("protein-containin...", "protein-cont. complex disassembly", df$Term)
  df$Term <- gsub("neuron apoptot...", "neuron apoptotic process", df$Term)
  df$Term <- gsub("benzene-containing compound metabolic pr...", "benzene-containing compound metabolic process", df$Term)
  df$Term <- gsub("multicellular organism gro...", "multicellular organism growth", df$Term)
  df$Term <- gsub("of circadian sleep/w...", "of circadian sleep/wake cycle, non-REM sleep", df$Term)
  df$Term <- gsub("of epithelial cell apical/ba...", "of epithelial cell apical/basal polarity", df$Term)
  df$Term <- gsub("acyl-CoA...", "acyl-CoA dehydrogenase", df$Term)
  df$Term <- gsub("proce...", "process", df$Term)
  df$Term <- gsub("prolif\\.\\.\\.", "proliferation", df$Term)
  df$Term <- gsub("cerebellar Purkinje cell-granule cell pr\\.\\.\\.",
                  "Purkinje cell-granule cell precursor signal", df$Term)
  df$Term <- gsub("smoothened signaling pathway involved in\\.\\.\\.",
                  "reg. of granule cell precursor cell prolif.", df$Term)
  df$Term <- gsub("cytokine-mediated\\.\\.\\.", "cytokine-med. sign. pathway", df$Term)
  df$Term <- gsub("metabolic pro\\.\\.\\.", "metabolism", df$Term)
  df$Term <- gsub("biosyntheti\\.\\.\\.", "biosynthesis", df$Term)
  df$Term <- gsub("oxidat\\.\\.\\.", "oxidation", df$Term)
  df$Term <- gsub("prolif\\.\\.\\.", "proliferation", df$Term)
  df$Term <- gsub("producti\\.\\.\\.", "production", df$Term)
  df$Term <- gsub("proliferat\\.\\.\\.", "proliferation", df$Term)

  df$Term <- gsub(" \\.\\.| \\.\\.\\.|, glu\\.\\.\\.| s\\.\\.\\.|pathw\\.\\.\\.", "", df$Term)
  df$Term <- gsub("humoral immune re...", "humoral immune response", df$Term)
  df$Term <- gsub("striated muscle c...", "striated muscle contraction", df$Term)
  df$Term <- gsub("endoplasmic retic...", "endoplasmic reticulum calcium ion concentration", df$Term)



  print(ggplot(df,
               aes(x = Significant/Annotated * 100, y = reorder(Term, -Fisher.elim_adjust),
                   color = Fisher.elim_adjust,
                   size = Significant))+
          geom_point()+
          ggtitle((add_title))+
          xlab("Hits (%)")+
          guides(color = guide_colorbar(order = 2, title = "Adj. p-value"),
                 size = guide_legend(order = 1, title = "DEG in GO"))+
          theme_minimal()+
          ylab("")+
          # theme(text = element_text(size = 28))  +
          theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1)) +
          theme(text = element_text(size = 26, color = "black"),
                title = element_text(size = 20),
                plot.title = element_text(hjust = 1),
                aspect.ratio = 1) +
          theme(axis.text = element_text(color = "black"))
  )
}


