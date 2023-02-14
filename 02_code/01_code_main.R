
# Code to generate all figures and analysis for the main part of the manuscript
# "Benefit of using interaction effects for the analysis of high-dimensional time-response or dose-response data for two-group comparisons"
#
# Contact: duda@statistik.tu-dortmund.de

# Setup ------------------------------------------------------------------------

# setwd("C:/Users/duda/Projekte/Git/gene_expression_interaction/02_Code")

# ONLY the required packages

library(ggplot2)
library(ggpubr)

source("./02_Code/00_functions.R")

# Data -------------------------------------------------------------------------

# Note:
# The data was provided by the IfADo and will be provided by request, given that the IfaDo allows it.


load("../01_data/010_gse.RData")

# Figure 1 ---------------------------------------------------------------------

df0 <- data.frame(x = c(1,2,1,2), y = c(0.2, 0.9, 0.7,1.4)* 10, col = factor(c("Group 0", "Group 0", "Group 1", "Group 1")))
p1 <- scheme_interaction_plot(df0, "a) interaction 0" )


df2 <- data.frame(x = c(1,2,1,2), y = c(0.2, 0.9, 0.7, 2)* 10, col = factor(c("Group 0", "Group 0", "Group 1", "Group 1")))
p2 <- scheme_interaction_plot(df2, "b) interaction +")


df <- data.frame(x = c(1,2,1,2), y = c(0.2, 0.9, 1.2,0.4)* 10, col = factor(c("Group 0", "Group 0", "Group 1", "Group 1")))
p3 <- scheme_interaction_plot(df, "c) interaction -")


df1 <- data.frame(x = c(1,2,1,2), y = c(0.2, 0.2, 0.2, 1.4)* 10, col = factor(c("Group 0", "Group 0", "Group 1", "Group 1")))
p4 <- scheme_interaction_plot(df1, "d) interaction +")


df3 <- data.frame(x = c(1,2,1,2), y = c(0.2, 0.9, 0.2, 1.5)* 10, col = factor(c("Group 0", "Group 0", "Group 1", "Group 1")))
p5 <- scheme_interaction_plot(df3, "e) interaction +")


df4 <- data.frame(x = c(1,2,1,2), y = c(0.2, 0.9, 2, 1.2) *10, col = factor(c("Group 0", "Group 0", "Group 1", "Group 1")))
p6 <- scheme_interaction_plot(df4, "f) interaction -")


# combine all plots
ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, common.legend = TRUE, legend="top")
ggsave(file = "./03_figures/01_main/fig_1_interaction_example_2.eps", device = cairo_ps, width = 10, height = 9)


# Figure 2, upper --------------------------------------------------------------

# dds object for method 1 

if(file.exists("dds_method1.RData")){ 
  load("020_dds_method1.RData")
}else{
  dds_method1 <- DESeqDataSet(gse_without_weeks, design = ~ diet + weeks)
  dds_method1 <- dds_method1[rowSums(counts(dds_method1)) >= 0.5, ]
  dds_method1 <- DESeq(dds_method1)
  #@caro : save einfuegen
}


plot_no_int("ENSMUSG00000069170", dds_method1, ref = 3, start = 3) +
  theme(legend.position = "none")
ggsave(file = "./03_figures/01_main/fig2_upper_example_gene.eps", width = 8, height = 4.5)
ggsave(file = "./03_figures/01_main/fig2_upper_example_gene.jpeg", width = 8, height = 4.5)

# Figure 2, lower --------------------------------------------------------------

# gse objekt for method ii, @caro namen ggf aendern? laden einfuegen?
gse_ia <- deseq_obj(gse_without_weeks, model = 2, threshold = 10, filter = T)


plot_it("ENSMUSG00000069170", gse_ia, ref = 3, start = 5.5, s1 = 1.5, s_arrow = 1.1, sl = 1.1)
ggsave(file = "./03_figures/01_main/fig2_lower_example_gene.eps", width = 8, height = 4.5)
ggsave(file = "./03_figures/01_main/fig2_lower_example_gene.jpeg", width = 8, height = 4.5)

# Figure 3 ---------------------------------------------------------------------

# create plots
pl_scheme_1 <- scheme_model(0.5, 0.9, 0.6, 1.11, 0.05, "not DEG", "not DEG", F, F, F, F,  1)
pl_scheme_2 <- scheme_model(0.5, 0.9, 0.6, 2.02, 0.05, "DEG", "DEG", F, T, F, T, 2)
pl_scheme_3 <- scheme_model(0.5, 1.3, 0.6, 2.5, 0.05, "not DEG", "DEG", T, T, T, T,  3)
pl_scheme_4 <- scheme_model(1.5, 2.47, 1.6, 0.71, 0.05, "not DEG", "DEG", T, T, T, T, 4)
pl_scheme_5 <- scheme_model(0.5, 1.01, 0.6, 1.43, 0.05, "DEG", "not DEG", F, T, F, F, 5)
pl_scheme_6 <- scheme_model(1.5, 1.89, 1.6, 0.55, 0.05, "DEG", "DEG", F, T, F, T, 6)
pl_scheme_7 <- scheme_model(1, 1.39, 1.1, 0.7, 0.05, "not DEG", "DEG", F, F, F, T, 7)


# save plots for method i and ii per scenario

plot_first_row <- ggarrange(pl_scheme_1[[1]],pl_scheme_1[[2]], ncol = 2, nrow = 1) 
annotate_figure(plot_first_row, fig.lab.size = 20,
                top = text_grob("Method I                   Method II", size = 25, face = "bold"))
ggsave(file = "./03_figures/01_main/fig3_scheme_1.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_1.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")

ggarrange(pl_scheme_2[[1]], pl_scheme_2[[2]],ncol = 2, nrow = 1)
ggsave(file = "./03_figures/01_main/fig3_scheme_2.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_2.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")

ggarrange(pl_scheme_3[[1]], pl_scheme_3[[2]],ncol = 2, nrow = 1)
ggsave(file = "./03_figures/01_main/fig3_scheme_3.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_3.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")

ggarrange(pl_scheme_4[[1]], pl_scheme_4[[2]],ncol = 2, nrow = 1)
ggsave(file = "./03_figures/01_main/fig3_scheme_4.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_4.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")

ggarrange(pl_scheme_5[[1]], pl_scheme_5[[2]],ncol = 2, nrow = 1)
ggsave(file = "./03_figures/01_main/fig3_scheme_5.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_5.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")

ggarrange(pl_scheme_6[[1]], pl_scheme_6[[2]],ncol = 2, nrow = 1)
ggsave(file = "./03_figures/01_main/fig3_scheme_6.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_6.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")

ggarrange(pl_scheme_7[[1]], pl_scheme_7[[2]],ncol = 2, nrow = 1, legend = "bottom", common.legend = T)
ggsave(file = "./03_figures/01_main/fig3_scheme_7.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")
ggsave(file = "./03_figures/01_main/fig3_scheme_7.jpeg", device = cairo_ps, width = 20, height = 10, unit  = "cm")


# Table 1 ----------------------------------------------------------------------
# results method i
if(file.exists(deg_mod_1_down)) #@caro: wie richtig? muss noch laden
else{
  # seperate dataset into the different weeks
  separate_dataset <- lapply(c(3,6), function(week){ #@caro: before: c(3,6,30, 36, 42, 48)
    # separate data set into each week
    gse_week       <- gse_without_weeks[,which(gse_without_weeks$weeks == week)]
    gse_week$weeks <- droplevels(gse_week$weeks)
    return(gse_week)}
  )
  
  # for each week, build model with ~ diet
  dds_obj_mod1    <- lapply(separate_dataset,
                            function(x) deseq_obj(gse = x,
                                                  model = 1,
                                                  threshold = 10, filter = T))
  
  # list of up and down regulated genes
  res_mod1        <- lapply(dds_obj_mod1,
                            function(x) de_mod1(dds = x, name = "diet_HFD_vs_SD"))
  # gene names for the degs
  gene_names_mod1 <- lapply(res_mod1, function(x) gene_names(x))
  
  # list separated into deg up and down reg genes
  mod_1_up <- lapply(1:length(gene_names_mod1),
                     function(x) c(gene_names_mod1[[x]]$up))
  mod_1_down <- lapply(1:length(gene_names_mod1),
                       function(x) c(gene_names_mod1[[x]]$down))
  # add names to all objects
  names(res_mod1) <- names(gene_names_mod1) <- names(mod_1_up) <-
    names(mod_1_down) <- names(dds_obj_mod1) <- paste0("week", c(3,6)) # @caro: before c(3,6,30, 36, 42, 48)
  
  # find "DEGs", in our case genes that are only DEG in week 6,
  # (p-value < 0.05 , abs(log2FC) > log2(1.5) in week 6;
  # p-value > 0.05 , abs(log2FC) < log2(1.5) in week 3), divided in up and down regulated
  deg_mod_1_down <- setdiff(mod_1_down$week6, c(mod_1_down$week3, mod_1_up$week3))
  deg_mod_1_up <- setdiff(mod_1_up$week6, c(mod_1_down$week3, mod_1_up$week3))
}


# table
# first row up, second down regulated
week_3_6_mod_1_table <- data.frame(
  "Week 3" = c(length(setdiff(mod_1_up$week3, c(mod_1_up$week6))), 
               length(setdiff(mod_1_down$week3, c(mod_1_down$week6)))),
  "Week 6" = c(length(setdiff(mod_1_up$week6, c(mod_1_up$week3))), 
               length(setdiff(mod_1_down$week6, c(mod_1_down$week3)))),
  "Overlap" = c(length(intersect(mod_1_up$week3, mod_1_up$week6)),
                length(intersect(mod_1_down$week3, mod_1_down$week6)))
)
week_3_6_mod_1_table

# Table 2 ----------------------------------------------------------------------
# first row up, second down regulated
data.frame("Method I only" = c(length(deg_mod_1_up) - sum(deg_mod_1_up %in% rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)]), 
                               c(length(deg_mod_1_down) - sum(deg_mod_1_down %in% rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)]))),
           "Overlap" = c(sum(deg_mod_1_up %in% rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)]),
                         sum(deg_mod_1_down %in% rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)])),
           "Method II only" = c(length(rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)]) - sum(rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)] %in% deg_mod_1_up),
                                length(rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)]) - sum(rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)] %in% deg_mod_1_down))
)
# Figure 4 ---------------------------------------------------------------------
plot_it("ENSMUSG00000025138", gse_ia, ref = 3, start = 9, s1 = 1.5, s_arrow = 1.1, sl = 1.1)

ggsave(file = "./03_figures/01_main/fig4_example_gene.eps", width = 8, height = 4.5)
ggsave(file = "./03_figures/01_main/fig4_example_gene.jpeg", width = 8, height = 4.5)
# Figure 5 ---------------------------------------------------------------------

#  method I: only the output of results = no shrinkage, no filtering @caro: warum hier no fitler und im naechsten shritt pre-filter
res_mod1_week3 <- results(dds_obj_mod1$week3, name = "diet_HFD_vs_SD")
res_mod1_week6 <- results(dds_obj_mod1$week6, name = "diet_HFD_vs_SD")


# because of pre-filter there are different genes in the objects, find intersection
intersect_names <- intersect(intersect(rownames(res_mod1_week3), rownames(res_mod1_week6)), rownames(res_mod2_ia))
res_mod1_week3_intersect <- data.frame(res_mod1_week3[intersect_names,])
res_mod1_week6_intersect <- data.frame(res_mod1_week6[intersect_names,])
res_mod2_ia_intersect <- data.frame(res_mod2_ia[intersect_names,])
res_mod2_main_intersect <- data.frame(res_mod2_main[intersect_names,])

# data frame with data from method ii for the points
df_no_sh <- data.frame("main" = res_mod2_main_intersect$log2FoldChange, 
                       "main_inter" = res_mod2_main_intersect$log2FoldChange + res_mod2_ia_intersect$log2FoldChange,  
                       "col" = 0)

#@caro ab hier weiter machen
# find "DEGs" for method i, in our case genes that are only DEG in week 6 (only by log2fc, not p-value)
deg_mod_1_no_sh_no_p   <-  which(abs(res_mod1_week6_intersect$log2FoldChange) > log2(1.5) &
                                   abs(res_mod1_week3_intersect$log2FoldChange) < log2(1.5))


# color for model 1: col = 1 for genes that are added in week 6 (up und down)
df_no_sh[deg_mod_1_no_sh_no_p, 3] <- 1

# interaktionseffekt sign (+2 -> die beides haben bekommen den wert 3)
# wert 2 ist nur method II, wert 3 ist method I and II, (res_mod2_ia is without shrinkage)
# +2 if interactino 
df_no_sh[which(abs(res_mod2_ia_intersect$log2FoldChange) > log2(1.5)), 3] <-
  (df_no_sh[which(abs(res_mod2_ia_intersect$log2FoldChange) > log2(1.5)), 3] + 2)


# change to factor and add names
df_no_sh$col <- as.factor(df_no_sh$col)
levels(df_no_sh$col) <- c("not significant", "Method I", "Method II", "Method I and II")

#split data set in sign and not sign (not sign points are plotted lighter)
df1_no_p_no_sh <- df_no_sh[df_no_sh$col != "not sign", ]
df2_no_p_no_sh <- df_no_sh[df_no_sh$col == "not sign", ]
ggplot(df1_no_p_no_sh, aes(x = main, y = main_inter, col = col)) + 
  geom_point(data = df2_no_p_no_sh, aes(x = main, y = main_inter), alpha = 0.05) +
  geom_point(alpha = 0.6) +
  #xlim(c(-3, 3)) +
  ylim(c(-3, 3)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()+ 
  ylab("main + interaction")+
  theme(axis.text.x = element_text(size = 22,angle = 45, hjust = 1))+
  theme(axis.text.y = element_text(size = 22))+
  theme(axis.title.x = element_text(size = 24)) + 
  theme(axis.title.y = element_text(size = 24))+ 
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 18))+ 
  #1
  annotate("text", x = 0.15, y = 0.51, label = 'bold("1")', parse = TRUE, size = 9) +
  geom_point(aes(x = 0.4, y = 0.51), col = "black") + 
  #2
  annotate("text", x = 0.15, y = 1.42, label = 'bold("2")', parse = TRUE, size = 9) +
  geom_point(aes(x = 0.4, y = 1.42), col = "black") + 
  #3
  annotate("text", x = 0.95, y = 1.9, label = 'bold("3")', parse = TRUE, size = 9) +
  geom_point(aes(x = 0.8, y = 1.9), col = "black") + 
  #4
  annotate("text", x = 1.15, y = -0.89, label = 'bold("4")', parse = TRUE, size = 9) +
  geom_point(aes(x = 0.97, y = -0.89), col = "black") + 
  #5
  annotate("text", x = 0.75, y = 0.83, label = 'bold("5")', parse = TRUE, size = 9) +
  geom_point(aes(x = 0.51, y = 0.83), col = "black") + 
  #6
  annotate("text", x = 0.15, y = -1.05, label = 'bold("6")', parse = TRUE, size = 9) +
  geom_point(aes(x = 0.39, y = -1.05), col = "black") + 
  #7
  annotate("text", x = 0.55,  y = -0.4, label = 'bold("7")', parse = TRUE, size = 9) + 
  geom_point(aes(x = 0.39, y = -0.4), col = "black") +
  geom_vline(xintercept = -log2(1.5), col = "black",  linetype = "dashed") +
  geom_vline(xintercept = log2(1.5), col = "black",  linetype = "dashed") +
  theme(legend.position="top") + 
  scale_x_continuous(breaks = c(round(-2, 0), -log2(1.5), 0, log2(1.5), 2),
                     labels = c(-2, "-log2(1.5)", 0, "log2(1.5)", 2), limits = c(-2,2))


ggsave(file = "./03_figures/01_main/fig5_scatter.eps", width = 25, height = 15, unit = "cm", 
       device = cairo_ps)
ggsave(file = "./03_figures/01_main/fig5_scatter.jpeg", width = 25, height = 15, unit = "cm", 
       device = cairo_ps)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# Suplementary Figure 1 --------------------------------------------------------

# Suplementary Figure 2 --------------------------------------------------------

# Suplementary Figure 3 --------------------------------------------------------











