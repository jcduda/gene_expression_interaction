
# Code to generate all figures and analysis for the main part of the manuscript
# "Benefit of using interaction effects for the analysis of high-dimensional time-response or dose-response data for two-group comparisons"
#
# Contact: duda@statistik.tu-dortmund.de

# Setup ------------------------------------------------------------------------

# setwd("C:/Users/duda/Projekte/Git/gene_expression_interaction/")

# ONLY the required packages

library(ggplot2)
library(ggpubr)

source("./02_Code/00_functions.R")

# Data -------------------------------------------------------------------------

# Note:
# The data was provided by the IfADo and will be provided by request, given that the IfaDo allows it.

load("./Data/010_gse.RData")

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


#ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, common.legend = TRUE, legend="bottom")
# combine all plots
ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, common.legend = TRUE, legend="top")
ggsave(file = "./03_figures/01_main/fig_1_interaction_example_2.eps", device = cairo_ps, width = 10, height = 9)
ggsave(file = "./03_figures/01_main/fig_1_interaction_example_2.jpeg", width = 10, height = 9)

# Figure 2, upper --------------------------------------------------------------

# ...

# Figure 2, lower --------------------------------------------------------------

# Figure 3 ---------------------------------------------------------------------

# Table 1 ----------------------------------------------------------------------

# Table 2 ----------------------------------------------------------------------

# Figure 4 ---------------------------------------------------------------------

# Figure 5 ---------------------------------------------------------------------



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------


# Suplementary Figure 1 --------------------------------------------------------

# Suplementary Figure 2 --------------------------------------------------------

# Suplementary Figure 3 --------------------------------------------------------











