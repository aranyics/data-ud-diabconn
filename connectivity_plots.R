# Required:
library(ggplot2)
library(reshape2)
library(scales)
library(cowplot)
library(RColorBrewer)
library(stringr)
# Optional:


group.connectivity = function(grp.data, grp.subjs, labels, grp1, grp2=NULL, lists=NULL)
{
  grp1.data = grp.data[,,which(sapply(grp.subjs, function(s) any(str_detect(s, grp1))))]
  if (!is.null(grp2))
    grp2.data = grp.data[,,which(sapply(grp.subjs, function(s) any(str_detect(s, grp2))))]
  
  
  grp1.mean = rowMeans(grp1.data, dims=2)
  grp1.sd = apply(grp1.data, c(1,2), sd)
  grp1.cv = grp1.sd / grp1.mean
  grp1.tstat = apply(grp1.data, c(1,2), function(m) t.test(m)$statistic)
  grp1.tstat.p = apply(grp1.data, c(1,2), function(m) t.test(m)$p.value)
  grp1.tstat.pfdr = array(p.adjust(grp1.tstat.p, 'fdr'), dim=c(length(labels), length(labels)))
  plot.grp1.mean = adjacency.plot(grp1.mean, labels, lists)
  plot.grp1.sd = adjacency.plot(grp1.sd, labels, lists)
  plot.grp1.cv = adjacency.plot(grp1.cv, labels, lists)
  plot.grp1.tstat = adjacency.plot(grp1.tstat, labels, lists)#, s.fdr=which(grp1.tstat.pfdr <= 0.05, arr.ind=TRUE))
  plot.grp1.str = strength.plot(grp1.mean, labels)
  
  if (!is.null(grp2))
  {
    grp2.mean = rowMeans(grp2.data, dims=2)
    grp2.sd = apply(grp2.data, c(1,2), sd)
    grp2.cv = grp2.sd / grp2.mean
    grp2.tstat = apply(grp2.data, c(1,2), function(m) t.test(m)$statistic)
    grp2.tstat.p = apply(grp2.data, c(1,2), function(m) t.test(m)$p.value)
    grp2.tstat.pfdr = array(p.adjust(grp2.tstat.p, 'fdr'), dim=c(length(labels), length(labels)))
    plot.grp2.mean = adjacency.plot(grp2.mean, labels, lists)
    plot.grp2.sd = adjacency.plot(grp2.sd, labels, lists)
    plot.grp2.cv = adjacency.plot(grp2.cv, labels, lists)
    plot.grp2.tstat = adjacency.plot(grp2.tstat, labels, lists)#, s.fdr=which(grp1.tstat.pfdr <= 0.05, arr.ind=TRUE))
    plot.grp2.str = strength.plot(grp2.mean, labels)
    
    grpmean = rowMeans(grp.data, dims=2)
    grpsd   = apply(grp.data, c(1,2), sd)
    grpcv   = grpsd / grpmean
    grpmeanratio = log(abs(grp1.mean / grp2.mean))
    grpsdratio = log(abs(grp1.sd / grp2.sd))
    grpcvratio = log(abs(grp1.cv / grp2.cv))
    grpmean.sym = grpmean + t(grpmean)
    grpmean.asym = grpmean - t(grpmean)
    grpdif.eucl = array(0, dim=c(length(labels), length(labels)))
    grpdif.tstat = array(0, dim=c(length(labels), length(labels)))
    grpdif.tstat.p = array(0, dim=c(length(labels), length(labels)))
    for (i in 1:length(labels))
      for (j in 1:length(labels))
      {
        grpdif.eucl[i, j] = sqrt(sum((grp1.data[i, j, ] - grp2.data[i, j, ]) ^ 2))
        grpdif.tstat[i, j] = t.test( grp1.data[i, j, ], grp2.data[i, j, ] )$statistic
        grpdif.tstat.p[i, j] = t.test( grp1.data[i, j, ], grp2.data[i, j, ] )$p.value
      }
    grpdif.tstat.pfdr = array(p.adjust(grpdif.tstat.p, 'fdr'), dim=c(length(labels), length(labels)))
    grpdif.tstat.pbonf = array(p.adjust(grpdif.tstat.p, 'bonferroni'), dim=c(length(labels), length(labels)))
    grpdif.tstat.padj = grpdif.tstat.pfdr
    
    p.unc = which(grpdif.tstat.p <= 0.001, arr.ind=TRUE)
    p.adj = which(grpdif.tstat.padj <= 0.05, arr.ind=TRUE)
    
    plot.grpmean = adjacency.plot(grpmean, labels, lists)
    plot.grpsd = adjacency.plot(grpsd, labels, lists)
    plot.grpcv = adjacency.plot(grpcv, labels, lists)
    plot.grpmeanratio = adjacency.plot(grpmeanratio, labels, lists)#, thresh.abs=log(2))
    plot.grpsdratio = adjacency.plot(grpsdratio, labels, lists)#, thresh.abs=log(2))
    plot.grpcvratio = adjacency.plot(grpcvratio, labels, lists)#, thresh.abs=log(2))
    plot.grpmean.sym = adjacency.plot(grpmean.sym, labels, lists)
    plot.grpmean.asym = adjacency.plot(grpmean.asym, labels, lists)
    plot.grpdif.eucl = adjacency.plot(grpdif.eucl, labels, lists)
    plot.grpdif.eucl.str = strength.plot(grpdif.eucl, labels)
    plot.grpdif.tstat = adjacency.plot(grpdif.tstat, labels, lists, s.unc=p.unc, s.fdr=p.adj)
  }
  
  if (is.null(grp2))
  {
    return ( list(plot.grp1.mean=plot.grp1.mean, plot.grp1.sd=plot.grp1.sd, plot.grp1.cv=plot.grp1.cv,
                  plot.grp1.tstat=plot.grp1.tstat,
                  plot.grp1.str=plot.grp1.str) )
  }
  else
  {  
    return ( list(plot.grp1.mean=plot.grp1.mean, plot.grp1.tstat=plot.grp1.tstat, plot.grp1.sd=plot.grp1.sd, plot.grp1.cv=plot.grp1.cv,
                  plot.grp2.mean=plot.grp2.mean, plot.grp2.tstat=plot.grp2.tstat, plot.grp2.sd=plot.grp2.sd, plot.grp2.cv=plot.grp2.cv,
                  plot.grp1.str=plot.grp1.str, plot.grp2.str=plot.grp2.str,
                  plot.grpmean=plot.grpmean, plot.grpmean.sym=plot.grpmean.sym, plot.grpmean.asym=plot.grpmean.asym,
                  plot.grpsd=plot.grpsd, plot.grpcv=plot.grpcv,
                  plot.grpmeanratio=plot.grpmeanratio,plot.grpsdratio=plot.grpsdratio, plot.grpcvratio=plot.grpcvratio,
                  plot.grpdif.tstat=plot.grpdif.tstat, plot.grpdif.eucl=plot.grpdif.eucl,
                  tstat.p.val=grpdif.tstat.p, tstat.p.val.adj=grpdif.tstat.padj,
                  plot.grpdif.eucl.str=plot.grpdif.eucl.str) )
  }
  
}


individual.connectivity = function(data, subjs, labels, lists, lim=NULL)
{
  plots = lapply(1:length(subjs), function(X){ adjacency.plot(data[,,X], labels, lists, lim=lim) })
  plots = setNames(plots, subjs)
  return(plots)
}


adjacency.plot = function(data, labels, lists=NULL, s.unc=NULL, s.fdr=NULL, lim=NULL, nw.lim=NULL, thresh.abs=NULL, nolegend=FALSE, puttext=FALSE)
{
  k = length(labels)
  
  #Full plot
  data.m = expand.grid( From=labels, To=labels )
  data.m$values = as.vector( t(data) )
  
  minval = round(min(data.m$values), 3)
  maxval = round(max(data.m$values), 3)
  #minval = round(quantile(c(data.m$values[which(-data.m$values > 0)]), probs=0.01), 3)
  #maxval = round(quantile(c(data.m$values[which(data.m$values > 0)]), probs=0.99), 3)
  #absmax = round(max(abs(data.m$values)), 2)
  if (!is.null(lim) && length(lim) == 2)
  {
    minval = lim[1]
    maxval = lim[2]
  }
  if (minval > 0 || is.na(minval)) minval = -0.001
  data.m$values[which(data.m$values < minval)] = minval
  data.m$values[which(data.m$values > maxval)] = maxval
  if (!is.null(thresh.abs))
  {
    data.m$values[which(abs(data.m$values) < thresh.abs)] = 0.0
  }
  
  g = ggplot(data.m, aes(From, To, fill = values)) + geom_tile() + coord_equal()
  if (puttext){
    g = g + geom_text(aes(label=round(values, 3)), colour="white", size=4.50)
    g = g + geom_text(aes(label=round(values, 3)), colour="white", size=4.51)
    g = g + geom_text(aes(label=round(values, 3)), colour="white", size=4.52)
    g = g + geom_text(aes(label=round(values, 3)), colour="white", size=4.52)
  }
  #g = g + scale_fill_gradient2(low = "blue", mid = "cyan", midpoint = 0, high = "red",
  #                             breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.01,maxval+0.01))
  g = g + scale_fill_gradientn(colours = c("#00f0ff", "#2415ff", "#2d004d", "#111111", "#990000", "#d28e00", "#ffff00"), #c("#D1E5F0", "#2166AC", "#111111", "#B2182B", "#FED976"), #rev(RColorBrewer::brewer.pal(7, "RdBu")),
                               values = rescale(c(minval-0.001, 0, maxval+0.001)),
                               breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.001,maxval+0.001))
  g = g + labs(fill=NULL)
  g = g + theme(panel.grid.major.x=element_blank(), #no gridlines
                panel.grid.minor.x=element_blank(), 
                panel.grid.major.y=element_blank(), 
                panel.grid.minor.y=element_blank(),
                panel.background=element_rect(fill="white"), # background=white
                axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 10, face = "bold"),
                plot.title = element_text(size = 20, face = "bold"),
                axis.text.y = element_text(size = 10,face = "bold"),
                legend.text = element_text(size = 10,face = "bold"))
  #legend.key.height = unit(0.5, units='in')  )
  g = g + scale_x_discrete(name = "", position = "top") + scale_y_discrete(name="", limits = rev(levels(data.m$To)))
  if (nolegend) g = g + theme(legend.position = "none")
  
  if (!is.null(lists))
  {
    for (l in unique(lists))
    {
      first = which(lists == l)[1]
      last = tail(which(lists == l),1)
      g = g + geom_rect(xmin = first-0.5, xmax = last+0.5, ymin = (k+1)-first+0.5, ymax = (k+1)-last-0.5, fill = NA, col = "white")
    }
    
    #Grouped plots - mean, pos.mean, neg.mean
    n = length(unique(lists))
    data2 = matrix(0, nrow=n, ncol=n)
    data2.pos = matrix(0, nrow=n, ncol=n)
    data2.neg = matrix(0, nrow=n, ncol=n)
    data[data==0] = NA
    for (lx in 1:n)
    {
      for (ly in 1:n)
      {
        data2[lx, ly] = mean( 
          as.matrix(data[ which(lists == unique(lists)[lx]), which(lists == unique(lists)[ly]) ]),
          na.rm=TRUE )
        if (is.nan(data2[lx, ly])) data2[lx, ly] = 0.0
        data2.pos[lx, ly] = mean( 
          data[ which(lists == unique(lists)[lx]), which(lists == unique(lists)[ly]) ][
            which(data[ which(lists == unique(lists)[lx]), which(lists == unique(lists)[ly]) ] > 0)
          ],
          na.rm=TRUE 
        )
        if (is.nan(data2.pos[lx, ly])) data2.pos[lx, ly] = 0.0
        data2.neg[lx, ly] = mean( 
          data[ which(lists == unique(lists)[lx]), which(lists == unique(lists)[ly]) ][
            which(data[ which(lists == unique(lists)[lx]), which(lists == unique(lists)[ly]) ] < 0)
          ],
          na.rm=TRUE 
        )
        if (is.nan(data2.neg[lx, ly])) data2.neg[lx, ly] = 0.0
      }
    }
    
    data2.m = expand.grid( From=unique(lists), To=unique(lists) )
    data2.m$values = as.vector( t(data2) )
    data2.pos.m = expand.grid( From=unique(lists), To=unique(lists) )
    data2.pos.m$values = as.vector( t(data2.pos) )
    data2.neg.m = expand.grid( From=unique(lists), To=unique(lists) )
    data2.neg.m$values = as.vector( t(data2.neg) )
    
    minval = round(min(data2.m$values), 4)
    maxval = round(max(data2.m$values), 4)
    if (!is.null(nw.lim) && length(nw.lim) == 2)
    {
      minval = nw.lim[1]
      maxval = nw.lim[2]
    }
    if (minval > 0 || is.na(minval)) minval = 0.0
    if (maxval < 0 || is.na(maxval)) maxval = 0.0
    data2.m$values[which(data2.m$values < minval)] = minval
    data2.m$values[which(data2.m$values > maxval)] = maxval
    
    gg = ggplot(data2.m, aes(From, To, fill = values)) + geom_tile() + coord_equal()
    if (puttext){
      gg = gg + geom_text(aes(label=round(values, 3)), colour="white", size=4.50)
      gg = gg + geom_text(aes(label=round(values, 3)), colour="white", size=4.51)
      gg = gg + geom_text(aes(label=round(values, 3)), colour="white", size=4.52)
      gg = gg + geom_text(aes(label=round(values, 3)), colour="white", size=4.52)
    }
    #gg = gg + scale_fill_gradient2(low = "red", mid = "black", midpoint = 0, high = "green",
    #                               breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.01,maxval+0.01))
    gg = gg + scale_fill_gradientn(colours = c("#00f0ff", "#2415ff", "#2d004d", "#111111", "#990000", "#d28e00", "#ffff00"),
                                   values = rescale(c(minval-0.001, 0, maxval+0.001)),
                                   breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.001,maxval+0.001))
    gg = gg + labs(fill=NULL)
    gg = gg + theme(panel.grid.major.x=element_blank(), #no gridlines
                    panel.grid.minor.x=element_blank(), 
                    panel.grid.major.y=element_blank(), 
                    panel.grid.minor.y=element_blank(),
                    panel.background=element_rect(fill="white"), # background=white
                    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 10, face = "bold"),
                    plot.title = element_text(size = 20, face = "bold"),
                    axis.text.y = element_text(size = 10,face = "bold"),
                    legend.text = element_text(size = 10,face = "bold"))
    #legend.key.height = unit(0.5, units='in')  )
    gg = gg + scale_x_discrete(name = "", position = "top") + scale_y_discrete(name="", limits = rev(levels(data2.m$To)))

    
    minval = 0.0
    maxval = round(max(data2.pos.m$values), 3)
    if (!is.null(nw.lim) && length(nw.lim) == 2)
    {
      maxval = nw.lim[2]
    }
    if (minval > 0 || is.na(minval)) minval = 0.0
    if (maxval < 0 || is.na(maxval)) maxval = 0.0
    data2.pos.m$values[which(data2.pos.m$values < minval)] = minval
    data2.pos.m$values[which(data2.pos.m$values > maxval)] = maxval
    
    ggp = ggplot(data2.pos.m, aes(From, To, fill = values)) + geom_tile() + coord_equal()
    #gg = gg + scale_fill_gradient2(low = "red", mid = "black", midpoint = 0, high = "green",
    #                               breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.01,maxval+0.01))
    ggp = ggp + scale_fill_gradientn(colours = c("#111111", "#990000", "#d28e00", "#ffff00"),
                                   values = rescale(c(0, maxval)),
                                   breaks = c(0,maxval), labels = c(0,maxval), limits = c(0.0, maxval))
    ggp = ggp + labs(fill=NULL)
    ggp = ggp + theme(panel.grid.major.x=element_blank(), #no gridlines
                    panel.grid.minor.x=element_blank(), 
                    panel.grid.major.y=element_blank(), 
                    panel.grid.minor.y=element_blank(),
                    panel.background=element_rect(fill="white"), # background=white
                    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 10, face = "bold"),
                    plot.title = element_text(size = 20, face = "bold"),
                    axis.text.y = element_text(size = 10,face = "bold"),
                    legend.text = element_text(size = 10,face = "bold"))
    #legend.key.height = unit(0.5, units='in')  )
    ggp = ggp + scale_x_discrete(name = "", position = "top") + scale_y_discrete(name="", limits = rev(levels(data2.pos.m$To)))
    
    
    minval = round(min(data2.neg.m$values), 3)
    maxval = 0.0
    if (!is.null(nw.lim) && length(nw.lim) == 2)
    {
      minval = nw.lim[1]
    }
    if (minval > 0 || is.na(minval)) minval = 0.0
    if (maxval < 0 || is.na(maxval)) maxval = 0.0
    data2.neg.m$values[which(data2.neg.m$values < minval)] = minval
    data2.neg.m$values[which(data2.neg.m$values > maxval)] = maxval
    
    ggn = ggplot(data2.neg.m, aes(From, To, fill = values)) + geom_tile() + coord_equal()
    #gg = gg + scale_fill_gradient2(low = "red", mid = "black", midpoint = 0, high = "green",
    #                               breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.01,maxval+0.01))
    ggn = ggn + scale_fill_gradientn(colours = c("#00f0ff", "#2415ff", "#2d004d", "#111111"),
                                   values = rescale(c(minval, 0)),
                                   breaks = c(minval,0), labels = c(minval,0), limits = c(minval, 0.0))
    ggn = ggn + labs(fill=NULL)
    ggn = ggn + theme(panel.grid.major.x=element_blank(), #no gridlines
                    panel.grid.minor.x=element_blank(), 
                    panel.grid.major.y=element_blank(), 
                    panel.grid.minor.y=element_blank(),
                    panel.background=element_rect(fill="white"), # background=white
                    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 10, face = "bold"),
                    plot.title = element_text(size = 20, face = "bold"),
                    axis.text.y = element_text(size = 10,face = "bold"),
                    legend.text = element_text(size = 10,face = "bold"))
    #legend.key.height = unit(0.5, units='in')  )
    ggn = ggn + scale_x_discrete(name = "", position = "top") + scale_y_discrete(name="", limits = rev(levels(data2.neg.m$To)))
    
    g = list(regions=g, networks=gg, networks.pos=ggp, networks.neg=ggn)
    
  }
  
  # Highlight significant
  if( !is.null(s.unc) && dim(s.unc)[1] > 0 )
  {
    for (i in 1:dim(s.unc)[1])
    {
      if (class(g) == 'list')
        g$regions = g$regions + geom_rect(xmin = s.unc[i,2]-0.5, xmax = s.unc[i,2]+0.5, ymin = (k+1)-s.unc[i,1]+0.5, ymax = (k+1)-s.unc[i,1]-0.5, fill = NA, col = "blue")
      else
        g = g + geom_rect(xmin = s.unc[i,2]-0.5, xmax = s.unc[i,2]+0.5, ymin = (k+1)-s.unc[i,1]+0.5, ymax = (k+1)-s.unc[i,1]-0.5, fill = NA, col = "blue")
    }
  }
  if( !is.null(s.fdr) && dim(s.fdr)[1] > 0 )
  {
    for (i in 1:dim(s.fdr)[1])
    {
      if (class(g) == 'list')
        g$regions = g$regions + geom_rect(xmin = s.fdr[i,2]-0.5, xmax = s.fdr[i,2]+0.5, ymin = (k+1)-s.fdr[i,1]+0.5, ymax = (k+1)-s.fdr[i,1]-0.5, fill = NA, col = "yellow")
      else
        g = g + geom_rect(xmin = s.fdr[i,2]-0.5, xmax = s.fdr[i,2]+0.5, ymin = (k+1)-s.fdr[i,1]+0.5, ymax = (k+1)-s.fdr[i,1]-0.5, fill = NA, col = "yellow")
    }
  }
  
  return(g)
}


bipolar.plot = function( data, lists, lim=NULL )
{
  library(dplyr)
  library(cowplot)
  
  pie.charts = list()
  
  n = length(unique(lists))
  for (lx in 1:n)
  {
    for (ly in 1:n)
    {      
      data.tmp = as.matrix(data[ which(lists == unique(lists)[lx]), which(lists == unique(lists)[ly]) ])
      
      minval = round(min(data.tmp), 2)
      maxval = round(max(data.tmp), 2)
      if (!is.null(lim) && length(lim) == 2)
      {
        minval = lim[1]#; data.tmp[which(data.tmp < minval)] = minval
        maxval = lim[2]#; data.tmp[which(data.tmp > maxval)] = maxval
      }
      if (minval > 0) minval = -0.01
      if (maxval < 0) maxval =  0.01
      
      chart.ID = str_c("C", lx, ly)
      plot.data = data.frame( values=c( length(which(data.tmp > 0)), length(which(data.tmp < 0)) ),
                              avg=c( mean(data.tmp[which(data.tmp > 0)]), mean(data.tmp[which(data.tmp < 0)]) ) )
      plot.data = plot.data %>%
        arrange(desc(values)) %>%
        mutate(prop = values / sum(plot.data$values) *100) %>%
        mutate(ypos = cumsum(prop) - 0.5*prop)
      plot.data = plot.data[which(plot.data$values != 0),]
      plot.data$pseud.avg = plot.data$avg
      plot.data$pseud.avg[which(plot.data$pseud.avg < minval)] = minval
      plot.data$pseud.avg[which(plot.data$pseud.avg > maxval)] = maxval
      
      g = ggplot(plot.data, aes(x="", y=prop, fill=pseud.avg)) + geom_bar(stat="identity", width=1, color="white")
      g = g + coord_polar("y", start=0)
      g = g + theme_void() + theme(legend.position="none")
      g = g + geom_text(aes(y = ypos, label = str_c(values, '\n', round(avg, 3))), color = "white", size=3)
      g = g + scale_fill_gradientn(colours = c("#00f0ff", "#2415ff", "#2d004d", "#111111", "#990000", "#d28e00", "#ffff00"),
                                   values = rescale(c(minval-0.01, 0, maxval+0.01)),
                                   breaks = c(minval,0,maxval), labels = c(minval,0,maxval), limits = c(minval-0.01,maxval+0.01))
      
      pie.charts[[ chart.ID ]] = g
    }
  }
  
  labels = as.list( c( unique(networks), rep('', (n-1)*n) ) )
  pie.merged = plot_grid( plotlist=pie.charts, nrow=n, labels=labels, label_size=8, scale=1.2 )
  
  return (pie.merged)
}


group.top.boxplot = function( grp.data, labels, grp.subjs, grp1, grp2, maxplots = 30, thresh.abs = NULL )
{
  grp1.data = grp.data[,,which(sapply(grp.subjs, function(s) any(str_detect(s, grp1))))]
  grp2.data = grp.data[,,which(sapply(grp.subjs, function(s) any(str_detect(s, grp2))))]
  
  grp1.mean = rowMeans(grp1.data, dims=2)
  grp1.sd   = apply(grp1.data, c(1,2), sd)
  grp1.cv   = grp1.sd / grp1.mean
  grp2.mean = rowMeans(grp2.data, dims=2)
  grp2.sd   = apply(grp2.data, c(1,2), sd)
  grp2.cv   = grp2.sd / grp2.mean
  #grpmean = rowMeans(grp.data, dims=2)
  #grpsd   = apply(grp.data, c(1,2), sd)
  
  grpmeanratio = log(grp1.mean / grp2.mean)
  grpsdratio = log(grp1.sd / grp2.sd)
  grpcvratio = log(grp1.cv / grp2.cv)
  
  if ( !is.null(thresh.abs) )
  {
    grpmeanratio[ which(abs(grpmeanratio) < thresh.abs) ] = 0.0
    grpsdratio[ which(abs(grpsdratio) < thresh.abs) ] = 0.0
    grpcvratio[ which(abs(grpcvratio) < thresh.abs) ] = 0.0
  }
  
  cat('mean ratio: min: ', min(grpmeanratio), ' max: ', max(grpmeanratio), ' threshold selected: ', length(which(grpmeanratio != 0.0)), '\n' )
  cat('sd ratio: min: ', min(grpsdratio), ' max: ', max(grpsdratio), ' threshold selected: ', length(which(grpsdratio != 0.0)), '\n' )
  cat('cv ratio: min: ', min(grpcvratio), ' max: ', max(grpsdratio), ' threshold selected: ', length(which(grpcvratio != 0.0)) )
  
  mean.idx = which(grpmeanratio != 0.0)
  sd.idx = which(grpsdratio != 0.0)
  cv.idx = which(grpcvratio != 0.0)
  if (length(mean.idx) > maxplots) mean.idx = mean.idx[order( abs(grpmeanratio[which(grpmeanratio != 0.0)]), decreasing=TRUE )][1:maxplots]
  if (length(sd.idx) > maxplots) sd.idx = sd.idx[order( abs(grpsdratio[which(grpsdratio != 0.0)]), decreasing=TRUE )][1:maxplots]
  if (length(cv.idx) > maxplots) cv.idx = cv.idx[order( abs(grpcvratio[which(grpcvratio != 0.0)]), decreasing=TRUE )][1:maxplots]
  
  mean.df = NULL
  sd.df = NULL
  cv.df = NULL
  
  # top sd diff
  for (i in 1:length(grp.subjs))
    for (j in 1:length(sd.idx))
    {
      id = grp.subjs[i]
      group = ''
      if (grp.subjs[i] %in% grp1 ) group = 'diabetes'
      if (grp.subjs[i] %in% grp2 ) group = 'obese'
      from = sd.idx[j] %/% length(labels) + 1
      to   = sd.idx[j] %%  length(labels)
      if ( (sd.idx[j] %% length(labels)) == 0 )
      {
        from = from - 1
        to = length(labels)
      }
      sd.df = rbind(sd.df, data.frame( id, group, str_c(labels[from], '-', labels[to]), grp.data[,,i][sd.idx[j]] ))
    }
  colnames(sd.df) = c('id', 'group', 'conn', 'value')
  
  # top cv diff
  for (i in 1:length(grp.subjs))
    for (j in 1:length(cv.idx))
    {
      id = grp.subjs[i]
      group = ''
      if (grp.subjs[i] %in% grp1 ) group = 'diabetes'
      if (grp.subjs[i] %in% grp2 ) group = 'obese'
      from = cv.idx[j] %/% length(labels) + 1
      to   = cv.idx[j] %%  length(labels)
      if ( (cv.idx[j] %% length(labels)) == 0 )
      {
        from = from - 1
        to = length(labels)
      }
      cv.df = rbind(cv.df, data.frame( id, group, str_c(labels[from], '-', labels[to]), grp.data[,,i][cv.idx[j]] ))
    }
  colnames(cv.df) = c('id', 'group', 'conn', 'value')
  
  # top mean diff
  for (i in 1:length(grp.subjs))
    for (j in 1:length(mean.idx))
    {
      id = grp.subjs[i]
      group = ''
      if (grp.subjs[i] %in% grp1 ) group = 'diabetes'
      if (grp.subjs[i] %in% grp2 ) group = 'obese'
      from = mean.idx[j] %/% length(labels) + 1
      to   = mean.idx[j] %%  length(labels)
      if ( (mean.idx[j] %% length(labels)) == 0 )
      {
        from = from - 1
        to = length(labels)
      }
      mean.df = rbind(mean.df, data.frame( id, group, str_c(labels[from], '-', labels[to]), grp.data[,,i][mean.idx[j]] ))
    }
  colnames(mean.df) = c('id', 'group', 'conn', 'value')
  
  # ggplots
  gg.sd = ggplot(sd.df, aes(x = conn, y = value, fill = group)) + geom_boxplot()
  gg.sd = gg.sd + geom_jitter(color="black", size=0.4, alpha=0.9)
  gg.cv = ggplot(cv.df, aes(x = conn, y = value, fill = group)) + geom_boxplot()
  gg.cv = gg.cv + geom_jitter(color="black", size=0.4, alpha=0.9)
  gg.mean = ggplot(mean.df, aes(x = conn, y = value, fill = group)) + geom_boxplot()
  gg.mean = gg.mean + geom_jitter(color="black", size=0.4, alpha=0.9)
  
  return( list(box.sd = gg.sd, box.cv = gg.cv, box.mean = gg.mean) )
}


relative.strength = function( data, labels, lists )
{
  data = abs(data)
  
  if (length(dim(data)) == 2)
  {
    data.i = array(0, dim=c(dim(data), 1))
    data.i[,,1] = data
  }
  else
    data.i = data
  rm(data)
  
  strength.out = matrix(0, nrow=dim(data.i)[3], ncol=length(labels), dimnames=list(NULL, labels))
  strength.in  = matrix(0, nrow=dim(data.i)[3], ncol=length(labels), dimnames=list(NULL, labels))
  for (s in 1:dim(data.i)[3])
  {
    for (i in 1:length(labels))
    {
      idx.in = which(lists == lists[i])
      idx.in = idx.in[-which(idx.in == i)]
      idx.ex = which(lists != lists[i])
      strength.out[s,i] = log( mean(data.i[idx.in,i,s]) / mean(data.i[idx.ex,i,s]) )
      strength.in[s,i]  = log( mean(data.i[i,idx.in,s]) / mean(data.i[i,idx.ex,s]) )
    }
  }
  
  strength = list(instr=strength.in, outstr=strength.out)
  
  return(strength)
}


strength.plot = function(data, label, lists=NULL)
{
  data = abs(data)
  
  strength = data.frame(label=label, In.strength=rowMeans(data), Out.strength=colMeans(data))
  if ( !is.null(lists) )
  {
    relstr = relative.strength(data, label, lists)
    strength = data.frame(label=label, In.strength=c(relstr$instr),
                          Out.strength=c(relstr$outstr))
  }
  
  plot.limits = c( min(min(strength$In.strength), min(strength$Out.strength)),
                   max(max(strength$In.strength), max(strength$Out.strength)) )
  plot.span = sqrt( (plot.limits[1]-plot.limits[2])^2 )
  plot.limits = c(plot.limits[1] - plot.span/10,
                  plot.limits[2] + plot.span/10 )
  
  line.slope = data.frame(sl = 1, int = 0, col = "green")
  
  g = ggplot(strength, aes(In.strength, Out.strength, label=label)) + geom_point(color=muted("red"), size=2) + coord_cartesian(xlim=plot.limits, ylim=plot.limits)
  g = g + geom_text(check_overlap=TRUE, nudge_y=plot.span/40, nudge_x=plot.span/25)
  g = g + geom_abline(data=line.slope, aes(slope=sl, intercept=int), size=1)
  g = g + theme_minimal()
  g = g + theme(axis.text.x = element_text(size = 12, face = "bold"),
                plot.title = element_text(size = 20, face = "bold"),
                axis.text.y = element_text(size = 12,face = "bold"),
                legend.text = element_text(size = 12,face = "bold"))
  g = g + geom_vline(xintercept = 0, col = "red")
  g = g + geom_hline(yintercept = 0, col = "red")
  
  return (g)
}


diverging.plot = function(data, labels)
{
  data[col(data) == row(data)] = 0
  
  regional.count = apply(array(1:length(labels)), 1, function(x){c(-1*(length(which(data[x,] < 0))+length(which(data[,x] < 0))),
                                                                   length(which(data[x,] > 0))+length(which(data[,x] > 0)))} )
  regional.weight = apply(array(1:length(labels)), 1, function(x){c(sum(c(data[x,which(data[x,] < 0)], data[which(data[,x] < 0),x])) / (2*length(labels)-1),
                                                                    sum(c(data[x,which(data[x,] > 0)], data[which(data[,x] > 0),x])) / (2*length(labels)-1))} )
  regional.weight[which(is.nan(regional.weight))] = 0.0
  
  lower.bound = quantile(regional.weight[1,], probs = 0.05)
  upper.bound = quantile(regional.weight[2,], probs = 0.95)
  
  df.count = data.frame(Region=labels, Inhibitory=regional.count[1,], Excitatory=regional.count[2,])
  df.count$Region = factor(df.count$Region, levels=rev(labels))
  df.count.m = melt(df.count)
  df.weight = data.frame(Region=labels, Inhibitory=regional.weight[1,], Excitatory=regional.weight[2,])
  df.weight$Region = factor(df.weight$Region, levels=rev(labels))
  df.weight.m = melt(df.weight)
  
  g1 = ggplot(df.count.m, aes(Region, value, fill=variable)) + geom_bar(stat='identity', width=0.5) + coord_flip()
  g1 = g1 + scale_fill_manual(values=c('#666666', '#999999')) + theme_bw() + theme(legend.position="none")
  g1 = g1 + ylab('Number of connections')
  
  g2 = ggplot(df.weight.m, aes(Region, value, fill=variable)) + geom_bar(stat='identity', width=0.5) + coord_flip()
  g2 = g2 + scale_fill_manual(values=c('#666666', '#999999')) + theme_bw() + theme(legend.position="none")
  g2 = g2 + ylab('Average connection strength difference')
  g2 = g2 + geom_hline(yintercept = lower.bound, color = "red", size = 1)
  g2 = g2 + geom_hline(yintercept = upper.bound, color = "red", size = 1)
  g2 = g2 + geom_hline(yintercept = 0, color = "black", size = 0.9)

  return ( list(count=g1, average=g2) )
}


flatten.connectivity = function(data, r.range, c.range, labelnames)
{
  data = data[r.range[1]:r.range[2], c.range[1]:c.range[2]]
  conn = array('', dim=(diff(r.range)+1) * (diff(c.range)+1))
  k = 1
  for (col in c.range[1]:c.range[2])
    for (row in r.range[1]:r.range[2])
    {
      conn[k] = str_c(labelnames[col], labelnames[row], sep=' --> ')
      k = k+1
    }
  df = data.frame(value=c(data), connection=conn, sign=sign(c(data)), stringsAsFactors=FALSE)
  df = df[order(df$value, decreasing=TRUE),]
  df$sign = factor(df$sign)
  return(df)
}


flatten.connectivity.idx = function(data, r.idx, c.idx, labelnames)
{
  data = data[r.idx, c.idx]
  conn = array('', dim=length(r.idx) * length(c.idx))
  k = 1
  for (col in c.idx)
    for (row in r.idx)
    {
      conn[k] = str_c(labelnames[col], labelnames[row], sep=' --> ')
      k = k+1
    }
  df = data.frame(value=c(data), connection=conn, sign=sign(c(data)), stringsAsFactors=FALSE)
  #df = df[order(df$value, decreasing=TRUE),]
  df$sign = factor(df$sign)
  return(df)
}


connection.barplot = function(data, labels, lists)
{
  labelnames = str_c(lists, labels, sep='_')
  plots = list()
  
  n = length(unique(lists))
  n.sizes = sapply(unique(lists), function(x) {length(which(lists == x))})
  n.pos = cumsum(n.sizes) - (n.sizes-1)
  for (lx in 1:n)
  {
    for (ly in 1:n)
    {
      label = str_c(unique(lists)[lx], unique(lists)[ly], sep='to')
      title = str_c('from ', unique(lists)[lx], ' to ', unique(lists)[ly])
      df = flatten.connectivity(data,
                                r.range=c(n.pos[ly], (n.pos[ly]+n.sizes[ly]-1)),
                                c.range=c(n.pos[lx], (n.pos[lx]+n.sizes[lx]-1)),
                                labelnames)
      g = ggplot(df, aes(reorder(connection, value), value, fill=sign)) + ggtitle(title)
      g = g + geom_bar(stat = 'identity') + coord_flip()
      g = g + geom_text(aes(connection, value, label=round(value, 4)), nudge_y=as.numeric(as.character(df$sign))*0.01, size=3)
      g = g + theme(panel.grid=element_blank(),
                    panel.background=element_rect(colour='white'),
                    axis.text.y=element_text(size=8, hjust=0),
                    legend.position='none')
      g = g + scale_y_discrete(name = NULL, limits=c(-0.05, 0.00, 0.05, 0.10, 0.15), expand=expand_scale(add=0.02))
      g = g + xlab(NULL)
      plots[[label]] = g
    }
  }
  
  return(plots)
}


connection.barplot.idx = function(data, labels, lists, threshold=-100.0)
{
  labelnames = str_c(lists, labels, sep='_')
  plots = list()
  
  n = length(unique(lists))
  n.sizes = sapply(unique(lists), function(x) {length(which(lists == x))})
  n.pos = cumsum(n.sizes) - (n.sizes-1)
  for (lx in 1:n)
  {
    label = str_c(unique(lists)[lx])
    title = str_c(unique(lists)[lx], ' endogenous')
    df = flatten.connectivity.idx(data,
                                  r.idx= n.pos[lx] : (n.pos[lx]+n.sizes[lx]-1),
                                  c.idx= n.pos[lx] : (n.pos[lx]+n.sizes[lx]-1),
                                  labelnames)
    df = df[which(df$value > threshold), ]
    g = ggplot(df, aes(reorder(connection, value), value, fill=sign)) + ggtitle(title)
    g = g + geom_bar(stat = 'identity') + coord_flip()
    g = g + geom_text(aes(connection, value, label=round(value, 4)), nudge_y=as.numeric(as.character(df$sign))*0.01, size=3)
    g = g + theme(panel.grid=element_blank(),
                  panel.background=element_rect(colour='white'),
                  axis.text.y=element_text(size=8, hjust=0),
                  legend.position='none')
    g = g + scale_y_discrete(name = NULL, limits=c(-0.05, 0.00, 0.05, 0.10, 0.15), expand=expand_scale(add=0.02))
    g = g + xlab(NULL)
    plots[[label]] = g
    
    label = str_c(unique(lists)[lx], '.out')
    title = str_c('from ', unique(lists)[lx])
    df = flatten.connectivity.idx(data,
                                  r.idx= seq(1:sum(n.sizes))[-c(n.pos[lx] : (n.pos[lx]+n.sizes[lx]-1))],
                                  c.idx= n.pos[lx] : (n.pos[lx]+n.sizes[lx]-1),
                                  labelnames)
    df = df[which(df$value > threshold), ]
    g = ggplot(df, aes(reorder(connection, value), value, fill=sign)) + ggtitle(title)
    g = g + geom_bar(stat = 'identity') + coord_flip()
    g = g + geom_text(aes(connection, value, label=round(value, 4)), nudge_y=as.numeric(as.character(df$sign))*0.01, size=3)
    g = g + theme(panel.grid=element_blank(),
                  panel.background=element_rect(colour='white'),
                  axis.text.y=element_text(size=8, hjust=0),
                  legend.position='none')
    g = g + scale_y_discrete(name = NULL, limits=c(-0.05, 0.00, 0.05, 0.10, 0.15), expand=expand_scale(add=0.02))
    g = g + xlab(NULL)
    plots[[label]] = g
    
    label = str_c(unique(lists)[lx], '.in')
    title = str_c('to ', unique(lists)[lx])
    df = flatten.connectivity.idx(data,
                                  r.idx= n.pos[lx] : (n.pos[lx]+n.sizes[lx]-1),
                                  c.idx= seq(1:sum(n.sizes))[-c(n.pos[lx] : (n.pos[lx]+n.sizes[lx]-1))],
                                  labelnames)
    df = df[which(df$value > threshold), ]
    g = ggplot(df, aes(reorder(connection, value), value, fill=sign)) + ggtitle(title)
    g = g + geom_bar(stat = 'identity') + coord_flip()
    g = g + geom_text(aes(connection, value, label=round(value, 4)), nudge_y=as.numeric(as.character(df$sign))*0.01, size=3)
    g = g + theme(panel.grid=element_blank(),
                  panel.background=element_rect(colour='white'),
                  axis.text.y=element_text(size=8, hjust=0),
                  legend.position='none')
    g = g + scale_y_discrete(name = NULL, limits=c(-0.05, 0.00, 0.05, 0.10, 0.15), expand=expand_scale(add=0.02))
    g = g + xlab(NULL)
    plots[[label]] = g
    
  }
  
  return(plots)
}


networks.masks = function(lists)
{
  sizes = sapply(unique(lists), function(x){length(which(lists == x))})
  idx = cumsum( sizes )
  idx = (idx - c(idx[1], diff(idx))) + 1
  template = array(0, dim=c(length(lists),length(lists)))
  masks = list(full = template + 1, empty = template)
  
  for (i in 1:length(sizes))
  {
    template = template * 0
    template[c(idx[i]:(idx[i]+sizes[i]-1)), c(idx[i]:(idx[i]+sizes[i]-1))] = 1
    masks[[str_c(names(sizes)[i], '.internal')]] = template
    template[c(idx[i]:(idx[i]+sizes[i]-1)), ] = 1
    template[ , c(idx[i]:(idx[i]+sizes[i]-1))] = 1
    masks[[str_c(names(sizes)[i], '.global')]] = template
  }
  
  return(masks)
}


regional.summary = function(data, labels)
{
  regstat = list()
  for (i in 1:length(labels))
  {
    regdata = c(data[i,-i], data[-i,i])
    regstat[[ labels[i] ]] = list( mean=mean(regdata))#, sd=sd(regdata) )
  }
  
  figdata = data.frame(regions=labels, avg=unlist(regstat))
  figdata$regions = factor(figdata$regions, levels=labels)
  
  g = ggplot(figdata, aes(x=regions, y=avg)) + geom_bar(stat='identity', width=0.8, color="#444444", fill="#CCCCCC")
  g = g + ylab("Average connectivity difference") + xlab(NULL)
  g = g + theme(panel.grid.major.x=element_blank(), #no gridlines
                panel.grid.minor.x=element_blank(), 
                panel.grid.major.y=element_blank(), 
                panel.grid.minor.y=element_blank(),
                panel.background=element_rect(fill="white"), # background=white
                axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0, size = 10, face = "bold"),
                plot.title = element_text(size = 20, face = "bold"),
                axis.text.y = element_text(size = 10,face = "bold"),
                legend.text = element_text(size = 10,face = "bold"))
  
  regstat[['figure']] = g
  
  return(regstat)
}


dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Dispersion test of count data:\n",
      length(x), " data points.\n",
      "Mean: ",mean(x),"\n",
      "Variance: ",var(x),"\n",
      "Probability of being drawn from Poisson distribution: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res)
}


dep.colnames = function( names )
{
  k = length(names)
  cnames = character()
  
  for ( i in 1:k )
  {
    for ( j in 1:k )
    {
      n = paste( names[i], names[j], sep='.' )
      
      cnames = c( cnames, n )
    }
  }
  
  return(cnames)
}
