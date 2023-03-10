
# Code to generate all figures and analysis for the main part of the manuscript
# "Benefit of using interaction effects for the analysis of high-dimensional time-response or dose-response data for two-group comparisons"
#
# Contact: duda@statistik.tu-dortmund.de

# Setup ------------------------------------------------------------------------

# setwd("C:/Users/duda/Projekte/Git/gene_expression_interaction/02_Code")

# ONLY the required packages

library(ggplot2)
library(ggpubr)
library(DESeq2)
library(MASS)

source("../02_Code/00_functions.R")

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
ggsave(file = "../03_figures/01_main/fig_1_interaction_example_2.eps", device = cairo_ps, width = 10, height = 9)


# Figure 2, upper --------------------------------------------------------------

# dds object for method 1 

if(file.exists("020_dds_method1.RData")){ 
  load("020_dds_method1.RData")
}else{
  dds_method1 <- DESeqDataSet(gse, design = ~ diet + weeks)
  dds_method1 <- dds_method1[rowSums(counts(dds_method1)) >= 0.5, ]
  dds_method1 <- DESeq(dds_method1)
  save(dds_method1, file = "020_dds_method1.RData")
}


plot_no_int("ENSMUSG00000069170", dds_method1, ref = 3, start = 3) +
  theme(legend.position = "none")
ggsave(file = "../03_figures/01_main/fig2_upper_example_gene.eps", width = 8, height = 4.5)
ggsave(file = "../03_figures/01_main/fig2_upper_example_gene.jpeg", width = 8, height = 4.5)

# Figure 2, lower --------------------------------------------------------------

# dds object method 2
if(file.exists("030_dds_method2.RData")){
  load("030_dds_method2.RData")
}else{
  dds_method2 <- deseq_obj(gse, model = 2, threshold = 10, filter = T)
  save(dds_method2, file = "030_dds_method2.RData")
}


plot_it("ENSMUSG00000069170", dds_method2, ref = 3, start = 5.5, s1 = 1.5, s_arrow = 1.1, sl = 1.1)
ggsave(file = "../03_figures/01_main/fig2_lower_example_gene.eps", width = 8, height = 4.5)
ggsave(file = "../03_figures/01_main/fig2_lower_example_gene.jpeg", width = 8, height = 4.5)

# Figure 3 ---------------------------------------------------------------------

# create plots for each scenario
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
ggsave(file = "../03_figures/01_main/fig3_scheme_1.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")


ggarrange(pl_scheme_2[[1]], pl_scheme_2[[2]],ncol = 2, nrow = 1)
ggsave(file = "../03_figures/01_main/fig3_scheme_2.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")


ggarrange(pl_scheme_3[[1]], pl_scheme_3[[2]],ncol = 2, nrow = 1)
ggsave(file = "../03_figures/01_main/fig3_scheme_3.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")


ggarrange(pl_scheme_4[[1]], pl_scheme_4[[2]],ncol = 2, nrow = 1)
ggsave(file = "../03_figures/01_main/fig3_scheme_4.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")


ggarrange(pl_scheme_5[[1]], pl_scheme_5[[2]],ncol = 2, nrow = 1)
ggsave(file = "../03_figures/01_main/fig3_scheme_5.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")


ggarrange(pl_scheme_6[[1]], pl_scheme_6[[2]],ncol = 2, nrow = 1)
ggsave(file = "../03_figures/01_main/fig3_scheme_6.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")


ggarrange(pl_scheme_7[[1]], pl_scheme_7[[2]],ncol = 2, nrow = 1, legend = "bottom", common.legend = T)
ggsave(file = "../03_figures/01_main/fig3_scheme_7.eps", device = cairo_ps, width = 20, height = 10, unit  = "cm")



# Table 1 ----------------------------------------------------------------------
# results method i
if(file.exists("040_res_mod1.RData")){
 load("040_res_mod1.RData") 
}else{
  # seperate dataset into the different weeks
  separate_dataset <- lapply(c(3,6), function(week){ 
    # separate data set into each week
    gse_week       <- gse[,which(gse$weeks == week)]
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
    names(mod_1_down) <- names(dds_obj_mod1) <- paste0("week", c(3,6)) 
  
  # find "DEGs", in our case genes that are only DEG in week 6,
  # (p-value < 0.05 , abs(log2FC) > log2(1.5) in week 6;
  # p-value > 0.05 , abs(log2FC) < log2(1.5) in week 3), divided in up and down regulated
  deg_mod_1_down <- setdiff(mod_1_down$week6, c(mod_1_down$week3, mod_1_up$week3))
  deg_mod_1_up <- setdiff(mod_1_up$week6, c(mod_1_down$week3, mod_1_up$week3))
  save(deg_mod_1_down, deg_mod_1_up, dds_obj_mod1, mod_1_down, mod_1_up, file = "040_res_mod1.RData")
}


# table
# first row up, second row down regulated
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
# find DEG with method ii
if(file.exists("050_results_model_2.RData")){
  load("050_results_model_2.RData")
}else{
  res_mod2 <- de_mod2("dietHFD.weeks6", dds_method2)
  save(res_mod2,  file = "050_results_model_2.RData")
}

# first row up, second down regulated
data.frame("Method I only" = c(length(deg_mod_1_up) - sum(deg_mod_1_up %in% rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)]), 
                               c(length(deg_mod_1_down) - sum(deg_mod_1_down %in% rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)]))),
           "Overlap" = c(sum(deg_mod_1_up %in% rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)]),
                         sum(deg_mod_1_down %in% rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)])),
           "Method II only" = c(length(rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)]) - sum(rownames(res_mod2)[which(res_mod2$log2FoldChange > 0)] %in% deg_mod_1_up),
                                length(rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)]) - sum(rownames(res_mod2)[which(res_mod2$log2FoldChange < 0)] %in% deg_mod_1_down))
)
# Figure 4 ---------------------------------------------------------------------
plot_it("ENSMUSG00000025138", dds_method2, ref = 3, start = 9, s1 = 1.5, s_arrow = 1.1, sl = 1.1)

ggsave(file = "../03_figures/01_main/fig4_example_gene.eps", width = 8, height = 4.5)
# Figure 5 ---------------------------------------------------------------------

#  method I: only the output of results = no shrinkage
res_mod1_week3 <- results(dds_obj_mod1$week3, name = "diet_HFD_vs_SD")
res_mod1_week6 <- results(dds_obj_mod1$week6, name = "diet_HFD_vs_SD")

# method II: get results without shrinkage for interaction effect at week 6 and main effect
res_mod2_ia <- results(dds_method2, name = "dietHFD.weeks6")
res_mod2_main <- results(dds_method2, name = "diet_HFD_vs_SD")


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


# find DEGs for method i, in our case genes that are only DEG in week 6 (only by log2fc, not p-value)
deg_mod_1_no_sh_no_p   <-  which(abs(res_mod1_week6_intersect$log2FoldChange) > log2(1.5) &
                                   abs(res_mod1_week3_intersect$log2FoldChange) < log2(1.5))


# color for model 1: col = 1 for genes that are added in week 6 (up and down)
df_no_sh[deg_mod_1_no_sh_no_p, 3] <- 1

# interaction effect significant: +2 -> col = 2:  method ii,
# col = 3: method i and ii
df_no_sh[which(abs(res_mod2_ia_intersect$log2FoldChange) > log2(1.5)), 3] <-
  (df_no_sh[which(abs(res_mod2_ia_intersect$log2FoldChange) > log2(1.5)), 3] + 2)


# change to factor and add names
df_no_sh$col <- as.factor(df_no_sh$col)
levels(df_no_sh$col) <- c("not significant", "Method I", "Method II", "Method I and II")

#split data set in sign and not sign 
df1_no_p_no_sh <- df_no_sh[df_no_sh$col != "not sign", ]
df2_no_p_no_sh <- df_no_sh[df_no_sh$col == "not sign", ]

ggplot(df1_no_p_no_sh, aes(x = main, y = main_inter, col = col)) + 
  geom_point(data = df2_no_p_no_sh, aes(x = main, y = main_inter), alpha = 0.05) +
  geom_point(alpha = 0.6) +
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


ggsave(file = "../03_figures/01_main/fig5_scatter.eps", width = 25, height = 15, unit = "cm", 
       device = cairo_ps)
ggsave(file = "../03_figures/01_main/fig5_scatter.jpeg", width = 25, height = 15, unit = "cm")
# Figure 6 ---------------------------------------------------------------------

# 1
#method i
p1.11 <- plot_it_3_6(g = "ENSMUSG00000037031", dds = dds_method2, ref = 3, start = 5.2, sl = 1.1, lim1 = 5.2, lim2 = 8.1, main = 0.1, interac = 0.1)+
  theme(legend.position="none")+
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank())
# method ii
p1.12 <- plot_it_3_6_mod1(g = "ENSMUSG00000037031", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 5.2,  sl = 1.1, lim1 = 5.2, lim2 = 8.1, main3 = 0.1, main6 = 0.1, number = paste("1: ", rowData(dds_method2["ENSMUSG00000037031",])[,2]))+
  theme(legend.position="none") +
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 


# 2
# method i
p2.21 <-  plot_it_3_6(g = "ENSMUSG00000035184", dds = dds_method2, ref = 3, start = 6, sl = 1.1, lim1 = 6, lim2 = 10.5, main = 0.08)+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 
# method ii
p2.22 <-   plot_it_3_6_mod1(g = "ENSMUSG00000035184", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 6,  sl = 1.1, lim1 = 6, lim2 = 10.5, main3 = 0.12, number = paste("2: ", rowData(dds_method2["ENSMUSG00000035184",])[,2]))+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 

# 3
# method i
p3.21 <- plot_it_3_6(g = "ENSMUSG00000022439", dds = dds_method2, ref = 3, start = 4,  sl = 1.1, lim1 = 4, lim2 = 8)+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 
# method ii
p3.22 <- plot_it_3_6_mod1(g = "ENSMUSG00000022439", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 4, sl = 1.1, lim1 = 4, lim2 = 8, number =paste("3: ", rowData(dds_method2["ENSMUSG00000022439",])[,2]))+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 

# 4
# method i
p4.21 <- plot_it_3_6(g = "ENSMUSG00000112013", dds = dds_method2, ref = 3, start = 2, sl = 1.1, lim1 = 2, lim2 = 6.5)+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 
# method ii
p4.22 <- plot_it_3_6_mod1(g = "ENSMUSG00000112013", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 2, s1 = 1.5,  sl = 1.1, lim1 = 2, lim2 = 6.5,number =paste("4: ", rowData(dds_method2["ENSMUSG00000112013",])[,2]))+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 

# 5
# method i
p5.11 <- plot_it_3_6(g = "ENSMUSG00000037348", dds = dds_method2, ref = 3, start = 10,  sl = 1.1, lim1 = 10, lim2 = 13, main = 0.1, interac = 0.1)+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 
# method ii
p5.12 <- plot_it_3_6_mod1(g = "ENSMUSG00000037348", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 10,  sl = 1.1, lim1 = 10, lim2 = 13, main3 =  0.1, number =paste("5: ", rowData(dds_method2["ENSMUSG00000037348",])[,2]))+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 

# 6
# method i
p6.11 <- plot_it_3_6(g = "ENSMUSG00000067951", dds = dds_method2, ref = 3, start = 4.5, sl = 1.1, lim1 = 4.5, lim2 = 7.5, main =  0.1)+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank())
# method ii
p6.12 <- plot_it_3_6_mod1(g = "ENSMUSG00000067951", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 4.5,  sl = 1.1, lim1 = 4.5, lim2 = 7.5, main3 = 0.1, number =paste("6: ", rowData(dds_method2["ENSMUSG00000067951",])[,2]))+  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 

# 7
# method i
p7.11 <- plot_it_3_6(g = "ENSMUSG00000028637", dds = dds_method2, ref = 3, start = 6.3,   sl = 1.1, lim1 = 6.3, lim2 = 8.2, main = 0.1)+ #+labs(x = "weeks", y = "")
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 
# method ii
p7.12 <- plot_it_3_6_mod1(g = "ENSMUSG00000028637", week3 = dds_obj_mod1[[1]], week6 = dds_obj_mod1[[2]], ref = 3, start = 6.3,  sl = 1.1, lim1 = 6.3, lim2 = 8.2, main3 = 0.1, main6 = -0.1, number =paste("7: ", rowData(dds_method2["ENSMUSG00000028637",])[,2]))+# +labs(x = "weeks", y = "")
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank()) 



# combining all plots 
plot_all <- ggarrange(p1.12, p1.11, p2.22, p2.21, p3.22, p3.21, p4.22, p4.21, 
                      p5.12, p5.11, p6.12, p6.11, p7.12, p7.11,
                      ncol = 2, nrow = 7, common.legend = TRUE, legend="bottom")

# add annotation
annotate_figure(plot_all, fig.lab.size = 20,
                top = text_grob("Method I             Method II", size = 25, face = "bold"),
                left = text_grob("Log2(normalized counts)", rot = 90, size = 25, face = "bold"))
ggsave("../03_figures/01_main/fig6_real_genes.eps", width = 15, height = 55, unit = "cm", 
       device = cairo_ps) 


# Supplementary Figure 1 --------------------------------------------------------

# object without filtering (50% threshold) the genes
if(file.exists("060_dds_method2_no_filter.Rdata")){
  load("060_dds_method2_no_filter.Rdata")
}else{
  dds_method2_no_filter <- deseq_obj(gse, model = 2, threshold = 10,
                                filter = F)
  save(dds_method2_no_filter, file = "060_dds_method2_no_filter.Rdata") 
}

# coefficient diet HDF vs SD
plot_hist(gse = dds_method2_no_filter, "diet_HFD_vs_SD")  + ggtitle("Coefficient WD vs SD")

ggsave(file = "../03_figures/02_supplement/fig1_hist_before_filter.eps", width = 5.5, height = 3.5)

# Supplementary Figure 2 --------------------------------------------------------

# coefficient diet HDF vs SD
plot_hist(gse = dds_method2, "diet_HFD_vs_SD") + ggtitle("Coefficient WD vs SD")

ggsave(file = "../03_figures/02_supplement/fig2_hist_after_filter.eps", width = 5.5, height = 3.5)

# Supplementary Figure 3 --------------------------------------------------------

# get shrinked results for main and interaction effect

if(file.exists("070_res_mod2_shrink.RData")){
  load("070_res_mod2_shrink.RData")
}else{
  res_mod2_ia_shrink <-  lfcShrink(dds_method2, coef = "dietHFD.weeks6", type = "apeglm")
  res_mod2_main_shrink <-  lfcShrink(dds_method2, coef = "diet_HFD_vs_SD", type = "apeglm")
  save(res_mod2_main_shrink, res_mod2_ia_shrink, file = "070_res_mod2_shrink.RData")
}

# create data.frame with coordinate values method ii with shrinkage

df_p <- data.frame("main" = res_mod2_main_shrink$log2FoldChange, 
                   "main_inter" = res_mod2_main_shrink$log2FoldChange +  res_mod2_ia_shrink$log2FoldChange, 
                   "significant" = abs(res_mod2_ia_shrink$log2FoldChange) > log2(1.5), 
                   "col" = 0)


# model 1 which have log foldchange < log(1.5 )

# color for model 1: col = 1 for genes that are added in week 6 (up und down)
df_p[deg_mod_1_down, 4] <- 1
df_p[deg_mod_1_up, 4] <- 1


# interaction effect significant: +2 -> col = 2:  method ii,
# col = 3: method i and ii
df_p[which(abs(res_mod2_ia_shrink$log2FoldChange) > log2(1.5) & res_mod2_ia_shrink$padj < 0.05), 4] <-
  (df_p[which(abs(res_mod2_ia_shrink$log2FoldChange) > log2(1.5)& res_mod2_ia_shrink$padj < 0.05), 4] + 2)

# change to factor
df_p$col <- as.factor(df_p$col)
levels(df_p$col) <- c("not significant", "Method I", "Method II", "Method I and II")

# divide into not sign and sign 
df1_p <- df_p[df_p$col != "not significant", ]
df2_p <- df_p[df_p$col == "not significant", ]

ggplot(df_p, aes(x = main, y = main_inter, col = col)) + 
  geom_point(alpha = 0.6) +
  theme_classic() +
  ylim(c(-3, 3)) +
  geom_abline(intercept = 0, slope = 1) +
  ylab("main + interaction")+
  theme(axis.text.x = element_text(size = 22))+
  theme(axis.text.y = element_text(size = 22))+
  theme(axis.title.x = element_text(size = 24)) + 
  theme(axis.title.y = element_text(size = 24))+ 
  theme(axis.text.x = element_text(size = 22,angle = 45, hjust = 1))+
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(legend.title = element_blank())+
  theme(legend.text = element_text(colour="black", size = 20)) +
  geom_vline(xintercept = -log2(1.5), col = "black",  linetype = "dashed") +
  geom_vline(xintercept = log2(1.5), col = "black",  linetype = "dashed") +
  
  #1
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000037031",2] - 0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000037031",2] + res_mod2_main_shrink["ENSMUSG00000037031",2], 
           label = 'bold("1")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000037031",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000037031",2] + res_mod2_main_shrink["ENSMUSG00000037031",2]), 
             col = "black") +  
  #2
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000035184",2] - 0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000035184",2] + res_mod2_main_shrink["ENSMUSG00000035184",2], 
           label = 'bold("2")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000035184",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000035184",2] + res_mod2_main_shrink["ENSMUSG00000035184",2]), 
             col = "black") + 
  #3
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000070323",2] + 0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000070323",2] + res_mod2_main_shrink["ENSMUSG00000070323",2], 
           label = 'bold("3")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000070323",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000070323",2] + res_mod2_main_shrink["ENSMUSG00000070323",2]), 
             col = "black") + 
  #4
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000112013",2] + 0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000112013",2] + res_mod2_main_shrink["ENSMUSG00000112013",2], 
           label = 'bold("4")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000112013",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000112013",2] + res_mod2_main_shrink["ENSMUSG00000112013",2]), 
             col = "black") + 
  #5
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000037348",2] +0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000037348",2] + res_mod2_main_shrink["ENSMUSG00000037348",2], 
           label = 'bold("5")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000037348",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000037348",2] + res_mod2_main_shrink["ENSMUSG00000037348",2]), 
             col = "black") + 
  #6
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000067951",2] - 0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000067951",2] + res_mod2_main_shrink["ENSMUSG00000067951",2], 
           label = 'bold("6")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000067951",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000067951",2] + res_mod2_main_shrink["ENSMUSG00000067951",2]), 
             col = "black") + 
  #7
  annotate("text", x = res_mod2_main_shrink["ENSMUSG00000028637",2] +0.2, 
           y = res_mod2_ia_shrink["ENSMUSG00000028637",2] + res_mod2_main_shrink["ENSMUSG00000028637",2], 
           label = 'bold("7")', parse = TRUE, size = 9) +
  geom_point(aes(x = res_mod2_main_shrink["ENSMUSG00000028637",2], 
                 y = res_mod2_ia_shrink["ENSMUSG00000028637",2] + res_mod2_main_shrink["ENSMUSG00000028637",2]), 
             col = "black") + 
  theme(legend.position="top") + 
  scale_x_continuous(breaks = c(-2, -log2(1.5), 0, 
                                log2(1.5), 2), 
                     labels = c(-2, "-log2(1.5)", 0,
                                "log2(1.5)", 2), 
                     limits = c(-2,2)) 


ggsave("../03_figures/02_supplement/fig3_scatter_annotated_W_p_shrink.eps", width = 25, height = 15, unit = "cm", 
       device = cairo_ps) 
ggsave("../03_figures/02_supplement/fig3_scatter_annotated_W_p_shrink.jpeg", width = 25, height = 15, unit = "cm") 







