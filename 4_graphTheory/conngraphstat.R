library(reshape2)
library(igraph)
library(stringr)
library(assortnet)
library(cowplot)


YOUR.SCRIPT.PATH = '/specify/your/script/path/4_graphtheory/'
setwd(YOUR.SCRIPT.PATH)

# Functions
source('connectivity_plots.R')
source('calc_stat.R')


# Files
D.csv = "data/BMR_1_diab_A.csv"
O.csv = "data/BMR_2_obes_A.csv"
DIFF.csv = "data/BMR_3_diff_A.csv"

D.ecm = as.matrix(read.csv(D.csv, header = FALSE))
O.ecm = as.matrix(read.csv(O.csv, header = FALSE))
DIFF.ecm = as.matrix(read.csv(DIFF.csv, header = FALSE))

DO.brain.csv = "results/stat_DO_brain.csv"
DO.nw.csv = "results/stat_DO_network.csv"
DO.reg.csv = "results/stat_DO_region.csv"
DIFF.brain.csv = "results/stat_DIFF_brain.csv"
DIFF.nw.csv = "results/stat_DIFF_network.csv"
DIFF.reg.csv = "results/stat_DIFF_region.csv"

write( str_c("matrix", "parameter", "Obesity", "T2DM", sep = ','), DO.brain.csv, append = FALSE)
write( str_c("matrix", "parameter", "interaction", "Obesity", "T2DM", sep=','), DO.nw.csv, append = FALSE)


# Region/network names
regions = c( as.matrix( read.csv("data/regionnames.csv", header = FALSE, stringsAsFactors = FALSE) ) )
networks = str_split_fixed(regions, '_', 2)[,1]
regions = str_split_fixed(regions, '_', 2)[,2]
nw = unique(networks)
RSN.comm = sapply(networks, function(x){ which(nw == x) })


# Reorder regions
regions_orig = c( as.matrix( read.csv("data/regionnames_correct.csv", header = FALSE, stringsAsFactors = FALSE) ) )
correct_order = sapply(regions_orig, function(x){which(regions == x)})

regions = regions[correct_order]
networks = networks[correct_order]

D.ecm = D.ecm[correct_order, correct_order]
O.ecm = O.ecm[correct_order, correct_order]
DIFF.ecm = DIFF.ecm[correct_order, correct_order]


# BMR plots
plot.obes = adjacency.plot(O.ecm, labels = regions, lists = networks, lim = c(-0.08, 0.32))
plot.diab = adjacency.plot(D.ecm, labels = regions, lists = networks, lim = c(-0.08, 0.32))
plot.diff = adjacency.plot(DIFF.ecm, labels = regions, lists = networks, lim = c(-0.05, 0.05))
plot_grid( plot.obes$regions, plot.diab$regions, plot.diff$regions,
           nrow = 1, labels = c('A - Obesity BMR', 'B - Diabetes BMR', 'C - Group difference BMR'))
ggsave('results/final_BMR.pdf', width = 21, height = 7)
plot.diff.diverge = diverging.plot(DIFF.ecm, regions)
plot_grid( plot.diff.diverge$average,
           nrow = 1, labels = NULL )
ggsave('results/final_BMR_diff_diverge.pdf', width = 8, height = 6)


# Connectivity matrix variations
O.diag = diag(O.ecm)
D.diag = diag(D.ecm)
DIFF.dig = diag(DIFF.ecm)

diag(O.ecm) = 0
diag(D.ecm) = 0
diag(DIFF.ecm) = 0

O.wd = abs(O.ecm)
D.wd = abs(D.ecm)
DIFF.wd = abs(DIFF.ecm)

O.wu = abs(O.ecm) + t(abs(O.ecm))
D.wu = abs(D.ecm) + t(abs(D.ecm))
DIFF.wu = abs(DIFF.ecm) + t(abs(DIFF.ecm))

O.sign = sign(O.ecm)
D.sign = sign(D.ecm)
DIFF.sign = sign(DIFF.ecm)


# REPRESENTATIVE MATRICES

# Diagonal
df = data.frame( reg = 1:length(O.diag), rsn = factor(RSN.comm), obes = O.diag, diab = D.diag )
M = melt(df, id.vars = c("reg", "rsn"))
write( str_c("ECM", "S(diag)", round(mean(O.diag),3), round(mean(D.diag),3), sep = ','),
       DO.brain.csv, append = TRUE)

df = data.frame( region=regions, T2DM=D.diag, Obesity=O.diag )
M = melt(df, id.vars = "region")
colnames(M) = c('Region', 'Group', 'S_diag')
M$Region = factor(M$Region,levels = regions)
ggplot(M, aes(x=Region, y=S_diag, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/region_ECM_Sdiag.pdf', width = 7, height = 4)

# WDCM - r
AO.wd = calc.stat(O.wd, RSN.comm, "assortativity")
AD.wd = calc.stat(D.wd, RSN.comm, "assortativity")
write( str_c("WDCM", "assortativity", round(AO.wd$r,3), round(AD.wd$r,3), sep = ','),
       DO.brain.csv, append = TRUE)

plot.r.O = adjacency.plot(AO.wd$mixing_matrix[1:7,1:7], nw, lim=c(0, 0.09), puttext = TRUE)
plot.r.D = adjacency.plot(AD.wd$mixing_matrix[1:7,1:7], nw, lim=c(0, 0.09), puttext = TRUE)
plot_grid(plot.r.D, plot.r.O, nrow = 1, labels = c('T2DM r', 'Obese r'))
ggsave('results/network_WDCM_r.pdf', width = 14, height = 7)

df = data.frame( network=nw, T2DM=AD.wd$nwout, Obesity=AO.wd$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/network7_WDCM_r.pdf', width = 3, height = 4)

# WDCM - Sin / Sout
gO.wd = graph_from_adjacency_matrix(O.wd, mode = "directed", weighted = TRUE)
gD.wd = graph_from_adjacency_matrix(D.wd, mode = "directed", weighted = TRUE)
SinO.wd = strength(gO.wd, mode = "in")
SinD.wd = strength(gD.wd, mode = "in")
SoutO.wd = strength(gO.wd, mode = "out")
SoutD.wd = strength(gD.wd, mode = "out")

df = data.frame( region=regions, T2DM=SinD.wd, Obesity=SinO.wd )
M = melt(df, id.vars = "region")
colnames(M) = c('Region', 'Group', 'S_in')
M$Region = factor(M$Region,levels = regions)
ggplot(M, aes(x=Region, y=S_in, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/region_WDCM_Sin.pdf', width = 7, height = 4)

df = data.frame( region=regions, T2DM=SoutD.wd, Obesity=SoutO.wd )
M = melt(df, id.vars = "region")
colnames(M) = c('Region', 'Group', 'S_out')
M$Region = factor(M$Region,levels = regions)
ggplot(M, aes(x=Region, y=S_out, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/region_WDCM_Sout.pdf', width = 7, height = 4)

# WDCM - balance
BO.wd = calc.stat(O.wd, RSN.comm, "balance")
BD.wd = calc.stat(D.wd, RSN.comm, "balance")
# write( str_c("WDCM", "balance", round(BO.wd$r,3), round(BD.wd$r,3), sep = ','),
#        DO.brain.csv, append = TRUE)
# write( str_c("WDCM", "balance(p)", round(BO.wd$rp,3), round(BD.wd$rp,3), sep = ','),
#        DO.brain.csv, append = TRUE)
write( str_c("WDCM", "balance -log10(p)", round(-log10(BO.wd$rp),3), round(-log10(BD.wd$rp),3), sep = ','),
       DO.brain.csv, append = TRUE)

plot.balp.O = adjacency.plot(BO.wd$outp, nw, lim=c(0,1), puttext = TRUE)
plot.balp.D = adjacency.plot(BD.wd$outp, nw, lim=c(0,1), puttext = TRUE)
BO.wd$outlogp = -log10(BO.wd$outp); BO.wd$outlogp[which(is.infinite(BO.wd$outlogp))] = 0
BD.wd$outlogp = -log10(BD.wd$outp); BD.wd$outlogp[which(is.infinite(BD.wd$outlogp))] = 0
plot.ballogp.O = adjacency.plot(BO.wd$outlogp, nw, lim=c(0, 5), puttext = TRUE)
plot.ballogp.D = adjacency.plot(BD.wd$outlogp, nw, lim=c(0, 5), puttext = TRUE)
plot.bal.O = adjacency.plot(BO.wd$out, nw, lim=c(-0.25, 0.5), puttext = TRUE)
plot.bal.D = adjacency.plot(BD.wd$out, nw, lim=c(-0.25, 0.5), puttext = TRUE)
plot_grid(plot.bal.D, plot.bal.O,
          plot.ballogp.D, plot.ballogp.O,
          nrow = 2, labels = c('T2DM balance', 'Obese balance',
                               '-log10(p)', '-log10(p)'))
ggsave('results/network_WDCM_balance.pdf', width = 14, height = 14)

df = data.frame( network=nw, T2DM=BD.wd$nwout, Obesity=BO.wd$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'Balance')
M$Network = factor(M$Network,levels = nw)
ggplot(M, aes(x=Network, y=Balance, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/network7_WDCM_balance.pdf', width = 3, height = 4)

df = data.frame( network=nw, T2DM=-log(BD.wd$nwoutp), Obesity=-log(BO.wd$nwoutp) )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'Balance_p')
M$Network = factor(M$Network,levels = nw)
ggplot(M, aes(x=Network, y=Balance_p, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + ylab("Balance (-log(p))") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/network7_WDCM_balance_-logp.pdf', width = 3, height = 4)

df = data.frame( region=regions, T2DM=c(BD.wd$regout), Obesity=c(BO.wd$regout) )
M = melt(df, id.vars = "region")
colnames(M) = c('Region', 'Group', 'Balance')
M$Region = factor(M$Region,levels = regions)
ggplot(M, aes(x=Region, y=Balance, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/region_WDCM_balance.pdf', width = 7, height = 4)

df = data.frame( region=regions, T2DM=c(-log(BD.wd$regoutp)), Obesity=c(-log(BO.wd$regoutp)) )
M = melt(df, id.vars = "region")
colnames(M) = c('Region', 'Group', 'Balance_p')
M$Region = factor(M$Region,levels = regions)
ggplot(M, aes(x=Region, y=Balance_p, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + ylab("Balance (-log(p))") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/region_WDCM_balance_-logpval.pdf', width = 7, height = 4)

df = data.frame( region=regions, T2DM=c(BD.wd$regoutp), Obesity=c(BO.wd$regoutp) )
M = melt(df, id.vars = "region")
colnames(M) = c('Region', 'Group', 'Balance_p')
M$Region = factor(M$Region,levels = regions)
ggplot(M, aes(x=Region, y=Balance_p, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + ylab("Balance (p)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/region_WDCM_balance_pval.pdf', width = 7, height = 4)

# WUCM/WDCM - strength
gO.wu = graph_from_adjacency_matrix(O.wu, mode = "undirected", weighted = TRUE)
gD.wu = graph_from_adjacency_matrix(D.wu, mode = "undirected", weighted = TRUE)
write( str_c("WUCM", "strength", round(mean(strength(gO.wu)),3), round(mean(strength(gD.wu)),3), sep = ','),
       DO.brain.csv, append = TRUE)

SO.wd = calc.stat(O.wd, RSN.comm, "strength")
SD.wd = calc.stat(D.wd, RSN.comm, "strength")
plot.str.O = adjacency.plot(SO.wd$mixing_matrix[1:7,1:7], nw, lim = c(0,8.068), puttext = TRUE)
plot.str.D = adjacency.plot(SD.wd$mixing_matrix[1:7,1:7], nw, lim = c(0,8.068), puttext = TRUE)
plot_grid(plot.str.D, plot.str.O, nrow = 1, labels = c('T2DM strength', 'Obese strength'))
ggsave('results/network_WDCM_strength.pdf', width = 14, height = 7)

df = data.frame( network=nw, T2DM=SD.wd$nwout, Obesity=SO.wd$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'Strength')
M$Network = factor(M$Network,levels = nw)
ggplot(M, aes(x=Network, y=Strength, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('results/network7_WDCM_strength.pdf', width = 3, height = 4)

# WUCM - SWP / dC / dL
SWP.diab.csv = read.csv('data/SWP_diab.csv', header=FALSE)
SWP.obes.csv = read.csv('data/SWP_obes.csv', header=FALSE)
write( str_c("WUCM", "SWP", mean(SWP.obes.csv$V1), mean(SWP.diab.csv$V1), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("WUCM", "dC", mean(SWP.obes.csv$V2), mean(SWP.diab.csv$V2), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("WUCM", "dL", mean(SWP.obes.csv$V3), mean(SWP.diab.csv$V3), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("WUCM", "SWP.sd", sd(SWP.obes.csv$V1), sd(SWP.diab.csv$V1), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("WUCM", "dC.sd", sd(SWP.obes.csv$V2), sd(SWP.diab.csv$V2), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("WUCM", "dL.sd", sd(SWP.obes.csv$V3), sd(SWP.diab.csv$V3), sep = ','),
       DO.brain.csv, append = TRUE)

# SDCM - D-/D+/S-/S+
write( str_c("SDCM", "D-", round(length(which(O.sign == -1)) / (36*36-36),3), round(length(which(D.sign == -1)) / (36*36-36),3), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("SDCM", "D+", round(length(which(O.sign == 1)) / (36*36-36),3), round(length(which(D.sign == 1)) / (36*36-36),3), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("SDCM", "S-", round(mean(O.wd[which(O.sign == -1)]),3), round(mean(D.wd[which(D.sign == -1)]),3), sep = ','),
       DO.brain.csv, append = TRUE)
write( str_c("SDCM", "S+", round(mean(O.wd[which(O.sign == 1)]),3), round(mean(D.wd[which(D.sign == 1)]),3), sep = ','),
       DO.brain.csv, append = TRUE)



# network masks
zeroconn = array(0, dim = c(length(regions), length(regions)))

diag.mask = zeroconn; diag.mask[row(diag.mask) == col(diag.mask)] = 1
intrinsic.masks = list()
extrinsic.masks = list()
total.masks = list()

global.intrinsic.mask = array(0, dim = c(length(regions), length(regions)))
global.extrinsic.mask = array(0, dim = c(length(regions), length(regions)))
global.total.mask = array(0, dim = c(length(regions), length(regions)))

for (i in 1:length(nw))
{
  nw.idx = which(networks == nw[i])
  empty.mask = zeroconn
  total.mask = empty.mask
  int.mask = empty.mask
  ext.mask = empty.mask
  total.mask[,nw.idx] = 1
  total.mask[nw.idx,] = 1
  int.mask[nw.idx, nw.idx] = 1
  ext.mask = total.mask - int.mask
  intrinsic.masks[[i]] = int.mask
  extrinsic.masks[[i]] = ext.mask
  total.masks[[i]] = total.mask
  global.intrinsic.mask = global.intrinsic.mask + int.mask
  global.extrinsic.mask = array(1 * as.logical(global.extrinsic.mask + ext.mask), dim = c(length(regions), length(regions)))
  global.total.mask = array(1 * as.logical(global.total.mask + total.mask), dim = c(length(regions), length(regions)))
}

rm(nw.idx, empty.mask, int.mask, ext.mask, total.mask)

diag(global.intrinsic.mask) = 0
O.wd.intra = O.wd * global.intrinsic.mask
O.wd.inter = O.wd * global.extrinsic.mask
D.wd.intra = D.wd * global.intrinsic.mask
D.wd.inter = D.wd * global.extrinsic.mask


# intra/inter network assortativity
AO.wd.intra = calc.stat(O.wd.intra, RSN.comm, "assortativity", sumgraph = sum(O.wd))
AD.wd.intra = calc.stat(D.wd.intra, RSN.comm, "assortativity", sumgraph = sum(D.wd))
df = data.frame( network=nw, Obesity=AO.wd.intra$nwout, T2DM=AD.wd.intra$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
gA.intra = ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Intra-network assortativity") +
  scale_fill_manual(values = c("#2415ff", "#d28e00"))
#ggsave('results/network7_WDCM_r_intra.pdf', width = 3, height = 4)
write( str_c("WDCM", "assortativity", "intra-network", round(mean(AO.wd.intra$nwout),3), round(mean(AD.wd.intra$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "assortativity", "intra-network", round(max(AO.wd.intra$nwout),3), round(max(AD.wd.intra$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "assortativity", "intra-network", round(min(AO.wd.intra$nwout),3), round(min(AD.wd.intra$nwout),3), sep=','), DO.nw.csv, append = TRUE)

AO.wd.inter = calc.stat(O.wd.inter, RSN.comm, "assortativity", sumgraph = sum(O.wd))
AD.wd.inter = calc.stat(D.wd.inter, RSN.comm, "assortativity", sumgraph = sum(D.wd))
df = data.frame( network=nw, Obesity=AO.wd.inter$nwout, T2DM=AD.wd.inter$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
gA.inter = ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Inter-network assortativity") +
  scale_fill_manual(values = c("#2415ff", "#d28e00"))
#ggsave('results/network7_WDCM_r_inter.pdf', width = 3, height = 4)
write( str_c("WDCM", "assortativity", "inter-network", round(mean(AO.wd.inter$nwout),3), round(mean(AD.wd.inter$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "assortativity", "inter-network", round(max(AO.wd.inter$nwout),3), round(max(AD.wd.inter$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "assortativity", "inter-network", round(min(AO.wd.inter$nwout),3), round(min(AD.wd.inter$nwout),3), sep=','), DO.nw.csv, append = TRUE)


# intra/inter network strength
SO.wd.intra = calc.stat(O.wd.intra, RSN.comm, "strength", sumgraph = sum(O.wd))
SD.wd.intra = calc.stat(D.wd.intra, RSN.comm, "strength", sumgraph = sum(D.wd))
df = data.frame( network=nw, Obesity=SO.wd.intra$nwout, T2DM=SD.wd.intra$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
gS.intra = ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Intra-network strength") +
  scale_fill_manual(values = c("#2415ff", "#d28e00"))
#ggsave('results/network7_WDCM_r_intra.pdf', width = 3, height = 4)
write( str_c("WDCM", "strength", "intra-network", round(mean(SO.wd.intra$nwout),3), round(mean(SD.wd.intra$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "strength", "intra-network", round(max(SO.wd.intra$nwout),3), round(max(SD.wd.intra$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "strength", "intra-network", round(min(SO.wd.intra$nwout),3), round(min(SD.wd.intra$nwout),3), sep=','), DO.nw.csv, append = TRUE)

SO.wd.inter = calc.stat(O.wd.inter, RSN.comm, "strength", sumgraph = sum(O.wd))
SD.wd.inter = calc.stat(D.wd.inter, RSN.comm, "strength", sumgraph = sum(D.wd))
df = data.frame( network=nw, Obesity=SO.wd.inter$nwout, T2DM=SD.wd.inter$nwout )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
gS.inter = ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Inter-network strength") +
  scale_fill_manual(values = c("#2415ff", "#d28e00"))
#ggsave('results/network7_WDCM_r_inter.pdf', width = 3, height = 4)
write( str_c("WDCM", "strength", "inter-network", round(mean(SO.wd.inter$nwout),3), round(mean(SD.wd.inter$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "strength", "inter-network", round(max(SO.wd.inter$nwout),3), round(max(SD.wd.inter$nwout),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "strength", "inter-network", round(min(SO.wd.inter$nwout),3), round(min(SD.wd.inter$nwout),3), sep=','), DO.nw.csv, append = TRUE)



# intra/inter network balance
BO.wd.intra = calc.stat(O.wd.intra, RSN.comm, "balance", sumgraph = sum(O.wd))
BD.wd.intra = calc.stat(D.wd.intra, RSN.comm, "balance", sumgraph = sum(D.wd))
df = data.frame( network=nw, Obesity=-log10(BO.wd.intra$nwoutp), T2DM=-log10(BD.wd.intra$nwoutp) )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
gB.intra = ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Intra-network balance") +
  scale_fill_manual(values = c("#2415ff", "#d28e00")) + ylim(c(0,8))
#ggsave('results/network7_WDCM_r_intra.pdf', width = 3, height = 4)
BO.wd.intra$nwoutp = -log10(BO.wd.intra$nwoutp)
BD.wd.intra$nwoutp = -log10(BD.wd.intra$nwoutp)
write( str_c("WDCM", "balance (-log10(p))", "intra-network", round(mean(BO.wd.intra$nwoutp),3), round(mean(BD.wd.intra$nwoutp),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "balance (-log10(p))", "intra-network", round(max(BO.wd.intra$nwoutp),3), round(max(BD.wd.intra$nwoutp),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "balance (-log10(p))", "intra-network", round(min(BO.wd.intra$nwoutp),3), round(min(BD.wd.intra$nwoutp),3), sep=','), DO.nw.csv, append = TRUE)

BO.wd.inter = calc.stat(O.wd.inter, RSN.comm, "balance", sumgraph = sum(O.wd))
BD.wd.inter = calc.stat(D.wd.inter, RSN.comm, "balance", sumgraph = sum(D.wd))
df = data.frame( network=nw, Obesity=-log10(BO.wd.inter$nwoutp), T2DM=-log10(BD.wd.inter$nwoutp) )
M = melt(df, id.vars = "network")
colnames(M) = c('Network', 'Group', 'r')
M$Network = factor(M$Network,levels = nw)
gB.inter = ggplot(M, aes(x=Network, y=r, fill=Group)) + geom_bar(stat = "identity", position = "dodge", width = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Inter-network balance") +
  scale_fill_manual(values = c("#2415ff", "#d28e00")) + ylim(c(0,8))
#ggsave('results/network7_WDCM_r_inter.pdf', width = 3, height = 4)
BO.wd.inter$nwoutp = -log10(BO.wd.inter$nwoutp)
BD.wd.inter$nwoutp = -log10(BD.wd.inter$nwoutp)
write( str_c("WDCM", "balance (-log10(p))", "inter-network", round(mean(BO.wd.inter$nwoutp),3), round(mean(BD.wd.inter$nwoutp),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "balance (-log10(p))", "inter-network", round(max(BO.wd.inter$nwoutp),3), round(max(BD.wd.inter$nwoutp),3), sep=','), DO.nw.csv, append = TRUE)
write( str_c("WDCM", "balance (-log10(p))", "inter-network", round(min(BO.wd.inter$nwoutp),3), round(min(BD.wd.inter$nwoutp),3), sep=','), DO.nw.csv, append = TRUE)



# Summary RSN figure
plot_grid( plot.str.O, plot.str.D,
           plot.r.O, plot.r.D,
           plot.ballogp.O, plot.ballogp.D,
           nrow=3, labels = c('A - Obesity Strength', 'B - T2DM Strength',
                              'C - Obesity Assortativity', 'D - T2DM Assortativity',
                              'E - Obesity Balance (-log10(p))', 'F - T2DM balance (-log10(p))'))
ggsave('results/final_network_summary1.pdf', width = 14, height = 21)

plot_grid( gS.intra, gS.inter, gA.intra, gA.inter, gB.intra, gB.inter,
           nrow = 3, labels = c('A', 'B', 'C', 'D', 'E', 'F'))
ggsave('results/final_network_summary2.pdf', width = 12, height = 15)

