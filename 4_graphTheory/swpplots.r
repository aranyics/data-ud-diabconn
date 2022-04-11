library(cowplot)
library(ggplot2)
library(reshape2)
library(stringr)

setwd('G:/PhD/diabetes/peb/')

diab1000 = read.csv('SWP_diab_1000.csv', header = FALSE)
obes1000 = read.csv('SWP_obes_1000.csv', header = FALSE)

swp = data.frame(T2DM=diab1000[,1], Obesity=obes1000[,1])
delta_C = data.frame(T2DM=diab1000[,2], Obesity=obes1000[,2])
delta_L = data.frame(T2DM=diab1000[,3], Obesity=obes1000[,3])

swp.m = melt(swp)
dC.m = melt(delta_C)
dL.m = melt(delta_L)

g1 = ggplot(swp.m, aes(x = variable, y = value, fill=variable)) + geom_boxplot() + ylab('SWP') + xlab(NULL)
g1 = g1 + theme(axis.text.x = element_text(size=12, face = "bold"),
                axis.title.y = element_text(size = 12, face = "bold"),
                legend.position = "none",
                plot.margin = unit(c(1, 1, 1, 1), "cm"))

g2 = ggplot(dC.m, aes(x = variable, y = value, fill=variable)) + geom_boxplot() + ylab('delta_C') + xlab(NULL)
g2 = g2 + theme(axis.text.x = element_text(size=12, face = "bold"),
                axis.title.y = element_text(size = 12, face = "bold"),
                legend.position = "none",
                plot.margin = unit(c(1, 1, 1, 1), "cm"))

g3 = ggplot(dL.m, aes(x = variable, y = value, fill=variable)) + geom_boxplot() + ylab('delta_L') + xlab(NULL)
g3 = g3 + theme(axis.text.x = element_text(size=12, face = "bold"),
                axis.title.y = element_text(size = 12, face = "bold"),
                legend.position = "none",
                plot.margin = unit(c(1, 1, 1, 1), "cm"))


plot_grid( g1, g2, g3,
           nrow = 1,
           labels = c('Small-world propensity', 'Clustering coefficient deviation', 'Path length deviation'))

ggsave('swp_1000.pdf', width = 16, height = 6)
