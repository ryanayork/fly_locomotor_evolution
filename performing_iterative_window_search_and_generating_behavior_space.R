############################
#####Initiate workspace#####
############################
rm(list=ls());
options(stringsAsFactors=F);
library(igraph);
source("TREBLE_functions.R")

#Set up colors for use in figures
dros_cols = c("coral", "red", "orangered3", "brown3", "darkred", 
              "red3", "indianred3", "salmon", "forestgreen",
              "burlywood", "deepskyblue2", "deepskyblue4", "gold2",
              rep('rosybrown3', 16))
names(dros_cols) = c("mauritiana", "santomea", "yakuba", "melanogaster", 
                     "simulans", "sechellia", "erecta", "teissieri",
                     "virilis", "willistoni", "pseudoobscura", "persimilis", "arizonae",
                     'la66', 'la69', 'z530', "z56", "z58", "zh16", "zh18", "zh20", "zh27",
                     "zh29", "zh32", "zh33", "zh34", "zh42", "zh47", "zh58")

###################
#####Load data#####
###################
#Load data
vel = readRDS('dros_locomotor_evolution_all_velocity_smoothed_xy_092520.RDS')

#Find which trials are longer than 1000 frames
vel = vel[which(unlist(lapply(vel, function(x) nrow(x)))>5000)]

#Sample 50 random trials
vel = vel[sample(seq(1, length(vel), 1), 100, replace = FALSE)]

#Simplify
vel = lapply(vel, function(x) data.frame(x = x$x,
                                         y = x$y,
                                         translational_velocity = x$vt_cm_smooth,
                                         angular_velocity = x$vr_cm_smooth,
                                         sideslip = x$vs_cm_smooth, 
                                         time = x$t,
                                         strain = x$strain,
                                         trial = x$exp,
                                         name = x$name))

#Combine
vel = do.call(rbind, vel)

#Split
vel = split(vel, vel$name)

#########################################
#####Perform iterative window search#####
#########################################
#Set up window sizes to test
toSweep = c(1, seq(2, 100, 2))

#Initialize list to save results from iterative_umap
iterative_windows = list()

#Run iterative windows and save results
for(i in 1:length(toSweep)){
  
  #Counter
  print(paste("Window size ", toSweep[i], '; ', i, ' out of ', length(toSweep), sep = ''))
  
  #Function
  iterative_windows[[as.character(toSweep[i])]] = iterative_umap(lapply(vel, function(x) x[1:2000,]),
                                                                 velocity_windows = TRUE,
                                                                 include_sideslip = TRUE, 
                                                                 plot = TRUE,
                                                                 filter_windows = TRUE,
                                                                 n_windows = 1000,
                                                                 max_ang_diff = 40,
                                                                 max_trans_diff = 0.04,
                                                                 max_side_diff = 0.04,
                                                                 symm = TRUE,
                                                                 window_size = toSweep[i])}

#Plot sample iterative umaps
par(mfrow = c(5,10), mar = c(1,1,1,1))
for(i in 1:length(iterative_windows)){
  plot(iterative_windows[[i]]$umaps[[40]][,1:2], 
       type = 'l')
}

#Save (if desired)
saveRDS(lapply(iterative_windows, function(x) x$umaps), 
        'fly_all_strains_50_individuals_iterative_windows_smooth_xy.RDS')

######################################################################
#####Analyze results of iterative window search (as in Figure S1)#####
######################################################################
#Remove umaps with less than n rows
for(i in 1:length(iterative_windows)){
  iterative_windows[[i]] = iterative_windows[[i]][lapply(iterative_windows[[i]], function(x) nrow(x))>999]
}

#Calculate Procrustes and Euclidean distance of results
iterative_windows_pr = lapply(iterative_windows, function(y) run_procrustes(y)$procrustes)
iterative_windows_dist = lapply(iterative_windows, function(y) run_procrustes(y)$euclidean_distances)

#Plot
par(mfrow = c(1,2))
plot_results(iterative_windows_pr,
             ylim = c(0, 10),
             ylab = "RMSE",
             xlab = "Window size (frames)",
             plot_as_lines = TRUE)
plot_results(iterative_windows_dist,
             ylim = c(0, 1000),
             ylab = "Mean Euclidean distance",
             xlab = "Window size (frames)",
             plot_as_lines = TRUE)

#Plot as variance
par(mfrow = c(1,2))
pr_var = plot_variance(iterative_windows_pr,
                       ylim = c(0, 0.2),
                       ylab = "RMSE",
                       xlab = "Window size (frames)",
                       return = TRUE)
euc_var = plot_variance(iterative_windows_dist,
                        ylim = c(0, 0.2),
                        ylab = "Mean Euclidean distance",
                        xlab = "Window size (frames)")

#As a line
par(mfrow = c(1,2))
plot(euc_var,
     xaxt = 'n',
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylab = "Euclidean distance (CV)",
     ylim = c(0,0.1),
     bty = 'n',
     las = 2,
     type = 'l',
     lwd = 2,
     col = 'grey50',
     xlab = "Window size (frames)")
axis(1, 
     at = seq(1, length(euc_var),1),
     labels = names(euc_var),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 2)

plot(pr_var,
     xaxt = 'n',
     cex.axis = 1.5,
     cex.lab = 1.5,
     ylab = "Procrustes distance (CV)",
     ylim = c(0,0.2),
     bty = 'n',
     las = 2,
     type = 'l',
     lwd = 2,
     col = 'grey50',
     xlab = "Window size (frames)")
axis(1, 
     at = seq(1, length(pr_var),1),
     labels = names(pr_var),
     cex.lab = 1.5,
     cex.axis = 1.5,
     las = 2)

#Calculate recurrence
recurrence = list()
for(i in 1:length(iterative_windows)){
  recurrence[[as.character(names(iterative_windows)[i])]] = calculate_recurrence(iterative_windows[[i]],
                                                                                 n_bins = 16)}

#Save (if desired)
saveRDS(recurrence, 
        'fly_all_strains_50_individuals_iterative_windows_recurrence.RDS')

#Plot recurrence distributions across replicates
plot_recurrence(recurrence,
                mar = c(0.25,0.25,0.25,0.25))

#Plot mean recurrence times by window size
means = lapply(recurrence, function(x) unlist(lapply(x, function(y) (median(y$recurrences)))))

par(mfrow = c(3,2))
plot_results(means,
             xlab = "Window size (frames)",
             ylab = "Mean recurrence time",
             ylim = c(0,200), 
             plot_as_lines = TRUE)

plot_variance(means,
              xlab = "Window size (frames)",
              ylab = "Coefficient of variation",
              ylim = c(0,1))

#Max bin
max_bin = lapply(recurrence, function(x) unlist(lapply(x, function(y) (which.max(y$histogram$counts)))))
plot_results(max_bin,
             ylim = c(0,200), 
             plot_as_lines = TRUE,
             ylab = "Maxiumum recurrent bin")

plot_variance(max_bin,
              xlab = "Window size (frames)",
              ylab = "Coefficient of variation",
              ylim = c(0,2))

#Percent recurrent
percent_recurrent = lapply(recurrence, function(x) unlist(lapply(x, function(y) y$total_proportion_recurrent)))
plot_results(percent_recurrent,
             ylim = c(0,1), 
             plot_as_lines = TRUE,
             ylab = "Total proportion recurrent")

plot_variance(percent_recurrent,
              xlab = "Window size (frames)",
              ylab = "Coefficient of variation",
              ylim = c(0,0.25))

#########################################
#####Extract windows of desired size#####
#########################################
#Load data
vel = readRDS('dros_locomotor_evolution_all_velocity_smoothed_xy_082820.RDS')

#Find which trials are longer than n frames
vel = vel[which(unlist(lapply(vel, function(x) nrow(x)))>300)]

#Interpolate to 30
time_res = 30
bandwidth = 0.25
for(i in 1:length(vel)){
  
  print(i)
  
  v = vel[[i]]
  
  #Interpolate
  len = max(v$t, na.rm = TRUE)-min(v$t, na.rm = TRUE)
  len_f = len*30
  step = len/len_f
  
  #Get intervals
  ints = seq(min(v$t, na.rm = TRUE), max(v$t, na.rm = TRUE), step)
  
  #Make df
  tmp = data.frame(x = approx(v$t, v$x, xout = ints, method = "linear")$y,
                   y = approx(v$t, v$y, xout = ints, method = "linear")$y,
                   angle = approx(v$t, v$angle, xout = ints, method = "linear")$y,
                   angle_r = approx(v$t, v$angle_r, xout = ints, method = "linear")$y,
                   theta = approx(v$t, v$theta, xout = ints, method = "linear")$y,
                   t = approx(v$t, v$t, xout = ints, method = "linear")$y)
  
  #Get diffs
  xdiff = diff(tmp$x)
  ydiff = diff(tmp$y)
  tdiff = diff(tmp$t)
  
  #Calculate vt
  print("Calculating vt")
  #tmp$vt = c(NA, (xdiff^2 + ydiff^2)*0.5)
  tmp$vt = c(NA, sqrt(xdiff^2 + ydiff^2))
  
  #Calculate vr (put over diff in time)
  print("Calculating vr")
  tmp$vr = c(NA, diff(tmp$theta))
  
  #Calculate vs
  #Equation from JAABA for sidweways velocity: (x mmt+1 − x mmt) cos(theta mmt + π/2) + (y mmt+1 − y mmt) sin(theta mmt + π/2)
  x2 = c(tmp$x[2:length(tmp$x)], NA)
  y2 = c(tmp$y[2:length(tmp$y)], NA)
  
  tmp$vs = (((x2-tmp$x)*cos(tmp$angle_r+(pi/2))) + ((y2-tmp$y)*sin(tmp$angle_r+(pi/2))))*100
  
  #Calculate acceleration (second order derivative of vt)
  print("Calculating acceleration")
  tmp$acc = c(NA, diff(tmp$vt))
  
  #Calculate vt and vr and vs at cm/s
  xdiff_cm = c(rep(0, time_res), diff(tmp$x, lag = time_res))
  ydiff_cm = c(rep(0, time_res), diff(tmp$y, lag = time_res))
  tdiff_cm = c(rep(0, time_res), diff(tmp$t, lag = time_res))
  adiff_cm = c(rep(0, time_res), diff(tmp$theta, lag = time_res))
  
  tmp$vt_cm = sqrt(xdiff_cm^2 + ydiff_cm^2)/tdiff_cm
  tmp$vr_cm = adiff_cm/tdiff_cm
  
  x2 = c(tmp$x[time_res:length(tmp$x)], rep(0, time_res-1))
  y2 = c(tmp$y[time_res:length(tmp$y)], rep(0, time_res-1))
  
  tmp$vs_cm = (((x2-tmp$x)*cos(tmp$angle_r+(pi/2))) + ((y2-tmp$y)*sin(tmp$angle_r+(pi/2))))
  tmp$acc_cm = c(rep(NA, time_res), diff(tmp$vt_cm, lag = time_res))
  
  tmp$vt_cm_smooth = tmp$vt_cm
  tmp$vr_cm_smooth = tmp$vr_cm
  tmp$vs_cm_smooth = tmp$vs_cm
  
  tmp$vt_cm_smooth[is.na(tmp$vt_cm_smooth)] = 0
  tmp$vr_cm_smooth[is.na(tmp$vr_cm_smooth)] = 0
  tmp$vs_cm_smooth[is.na(tmp$vs_cm_smooth)] = 0
  
  tmp$vt_cm_smooth = ksmooth(tmp$t, tmp$vt_cm_smooth, kernel = "normal", bandwidth = bandwidth)$y
  tmp$vr_cm_smooth = ksmooth(tmp$t, tmp$vr_cm_smooth, kernel = "normal", bandwidth = bandwidth)$y
  tmp$vs_cm_smooth = ksmooth(tmp$t, tmp$vs_cm_smooth, kernel = "normal", bandwidth = bandwidth)$y
  
  #Simplify
  vel[[i]] = data.frame(x = tmp$x,
                        y = tmp$y,
                        translational_velocity = tmp$vt_cm_smooth,
                        angular_velocity = tmp$vr_cm_smooth,
                        sideslip = tmp$vs_cm_smooth, 
                        time = tmp$t,
                        strain = rep(v$strain[1], nrow(tmp)),
                        trial = rep(v$exp[1], nrow(tmp)),
                        name = rep(v$name[1], nrow(tmp)))  
}

#Get windows
win = list()
for(i in 1:length(vel)){
  print(paste(i, 'out of', length(vel)))
  
  win[[names(vel)[i]]] = get_velocity_windows(vel[[i]],
                                              include_sideslip = TRUE, 
                                              filter_windows = TRUE,
                                              max_ang_diff = 40,
                                              max_trans_diff = 0.1,
                                              max_side_diff = 0.1,
                                              symm = TRUE,
                                              name = names(vel)[i],
                                              window_size = 10)
  
  saveRDS(win[[names(vel)[i]]], 
          paste('../../01_analysis_files/00_windows/',
                names(vel[i]),
                '_velocity_windows_092220_30hz_size10.RDS',
                sep = ''))}

#Combine
wins = do.call(cbind, win)

#Save (if desired)
saveRDS(wins, 'all_trials_velocity_windows_30hz_size10.RDS')

##########################################
#####Generate behavior space via UMAP#####
##########################################
#Run UMAP
u = umap(t(wins), verbose = TRUE)

#Get layout
layout = data.frame(x = u$layout[,1],
                    y = u$layout[,2])

#Save
saveRDS(layout, 'all_trials_umap_layout_30hz_size10.RDS')
saveRDS(u, 'all_trials_umap_30hz_size10.RDS')

##############################################################
#####Characterize behavior space (as in Figures 2 and S1)#####
##############################################################
#Load layout
layout = readRDS("all_trials_umap_layout_30hz_size10.RDS")

#Load windows
win = readRDS('all_trials_velocity_windows_30hz_size10.RDS')

#Bin
layout = bin_umap(layout, n_bins = 64)$layout

#Plot behavior space
par(mfrow = c(1,2))
plot(layout[,1:2],
     pch = 20,
     col = alpha('gray50', 0.05),
     cex = 0.05,
     bty = 'n',
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '')
plot(layout[,1:2],
     col = alpha('gray50', 0.05),
     type = 'l',
     bty = 'n',
     lwd = 0.1,
     xaxt = 'n',
     yaxt = 'n',
     xlab = '',
     ylab = '')
dev.off()

#Plot vector field
vector_field = plot_vector_field(layout, return = TRUE)

#Plot colored by vector lengths
d = round(vector_field$dist, 2)
cols = rev(viridis_pal(option = 'A')(length(seq(0, max(d), 0.01))))
names(cols) = seq(0, max(d), 0.01)
cols = cols[match(d, names(cols))]

par(mar = c(1,1,1,1))
plot(unlist(lapply(strsplit(rownames(vector_field), "_"), function(v){v[1]})),
     unlist(lapply(strsplit(rownames(vector_field), "_"), function(v){v[2]})),
     pch = 20,
     col = cols,
     bty = 'n',
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '')

#Plot velocity features
plot_umap_features(layout, 
                   win, 
                   n_features = 3,
                   bin_umap = TRUE,
                   n_bins = 200,
                   cex = 0.25,
                   colors = c('darkgreen', 'darkmagenta', 'gold4'))

#Add time
layout$time = lapply(strsplit(colnames(win), "_"), function(v){v[1]})

#Add strain
layout$strain = lapply(strsplit(colnames(win), "_"), function(v){v[3]})

#Add individual
layout$individual = paste(lapply(strsplit(colnames(win), "_"), function(v){v[3]}),
                          lapply(strsplit(colnames(win), "_"), function(v){v[4]}),
                          lapply(strsplit(colnames(win), "_"), function(v){v[5]}),
                          lapply(strsplit(colnames(win), "_"), function(v){v[6]}),
                          lapply(strsplit(colnames(win), "_"), function(v){v[7]}),
                          sep = '_')

#Save
saveRDS(layout, 
        "all_trials_umap_layout_annotated_30hz_size10.RDS")

#Plot n trials per species
inds = unique(layout$individual)
inds = table(unlist(lapply(strsplit(inds, "_"), function(v){v[1]})))
inds = sort(inds)
cols = dros_cols[match(names(inds), names(dros_cols))]
plot(inds, 
     pch = 20, 
     col = cols,
     las = 2,
     ylab = 'n trials')

##############################################
#####Louvain clustering (as in Figure S2)#####
##############################################
#Rebin
layout = bin_umap(layout, n_bins = 64)$layout

#Make into graph
g = graph_from_data_frame(as.data.frame(
  cbind(layout[1:(nrow(layout)-1),]$xy_new, 
        layout[2:nrow(layout),]$xy_new)), directed = FALSE)

#Cluster with Louvain
partition <- cluster_louvain(g)

#Plot membership
m = partition$membership
names(m) = partition$names
cols = WGCNA::standardColors(length(partition))[m]

plot(unlist(lapply(strsplit(names(m), "_"), function(v){v[1]})),
     unlist(lapply(strsplit(names(m), "_"), function(v){v[2]})),
     pch = 20,
     col = cols,
     bty = 'n',
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '')

#Add to layout
layout$louvain_cluster = rep(NA, nrow(layout))
for(i in 1:length(m)){
  print(i)
  layout[layout$xy_new == names(m)[i],]$louvain_cluster = m[i]
}

##Mean velocities per cluster
#Load windows
win = readRDS('all_trials_velocity_windows_30hz_size10.RDS')

#Add velocities to layout
layout$vt = colMeans(win[1:10,])
layout$vr = colMeans(win[11:20,])
layout$vs = colMeans(win[21:30,])

#Save
saveRDS(layout,
        'all_trials_umap_layout_annotated_30hz_size10_with_louvain_clusters.RDS')

#Ternary plot
cols = RColorBrewer::brewer.pal(length(unique(layout$louvain_cluster)), 'Paired')
names(cols) = unique(layout$louvain_cluster)

library(Ternary)
dat = layout[,11:13]
dat = apply(dat, 2, function(x) x/max(x))

for(i in 1:7){
  d = dat[layout$louvain_cluster == i,]
  
  png(paste('cluster', i, '_ternary_plot.png', sep = ''), width = 800, height = 800)
  TernaryPlot(axis.labels = seq(0,0,0), 
              alab = '', 
              blab = '', 
              clab = '',
              grid.lines=0, 
              grid.lty='dotted',
              grid.minor.lines=1,
              grid.minor.lty='dotted',
              lab.cex = 1.5,
              axis.cex = 1.5)
  
  TernaryPoints(d, col = alpha(cols[i], 0.25), pch = '.')
  segments(-0.5, 0, 0, 0.865, lwd = 6)
  segments(0.5, 0, 0, 0.865, lwd = 6)
  segments(-0.5, 0, 0.5, 0, lwd = 6)
  
  dev.off()
}

#Split
cl_vels = split(layout, layout$louvain_cluster)

#Calculate mean and variance
vt_m = unlist(lapply(cl_vels, function(x) median(x$vt)))
vr_m = unlist(lapply(cl_vels, function(x) median(x$vr)))
vs_m = unlist(lapply(cl_vels, function(x) median(x$vs)))

vt_e = lapply(cl_vels, function(x) boxplot.stats(x$vt)$stats)
vr_e = lapply(cl_vels, function(x) boxplot.stats(x$vr)$stats)
vs_e = lapply(cl_vels, function(x) boxplot.stats(x$vs)$stats)

#Plot
plot(vt_m,
     seq(1, length(vt_m), 1),
     type = 'n',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlim = c(0,1.5),
     bty = 'n',
     yaxt = 'n',
     xlab = 'Translational velocity (cm/sec)',
     ylab = '')
segments(unlist(lapply(vt_e, function(x) x[2])),
         seq(1, length(vt_m), 1),
         unlist(lapply(vt_e, function(x) x[4])),
         seq(1, length(vt_m), 1),
         col = cols,
         lwd = 1.5)
points(vt_m,
       seq(1, length(vt_m), 1),
       pch = 20,
       col = cols,
       cex = 3)

plot(vr_m,
     seq(1, length(vr_m), 1),
     type = 'n',
     cex.axis = 1.5,
     cex.lab = 1.5,
     xlim = c(0,100),
     bty = 'n',
     xlab = 'Angular velocity (degrees/second)',
     ylab = 'Cluster')
segments(unlist(lapply(vr_e, function(x) x[2])),
         seq(1, length(vr_m), 1),
         unlist(lapply(vr_e, function(x) x[4])),
         seq(1, length(vr_m), 1),
         col = cols,
         lwd = 1.5)
points(vr_m,
       seq(1, length(vr_m), 1),
       pch = 20,
       col = cols,
       cex = 3)

#Plot as individual points
partition = unique(layout$louvain_cluster)
m = layout$louvain_cluster[1:20000]
cols = RColorBrewer::brewer.pal(length(unique(layout$louvain_cluster)), 'Paired')[m]
plot(layout[1:20000,1],
     layout[1:20000,2],
     pch = 20,
     col = alpha(cols, 0.5),
     cex = 0.25,
     bty = 'n',
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '')

#Calculate centers
centers = list()
for(i in 1:length(partition)){
  x = layout$x[layout$louvain_cluster == unique(layout$louvain_cluster)[i]]
  y = layout$y[layout$louvain_cluster == unique(layout$louvain_cluster)[i]]
  centers[[as.character(unique(layout$louvain_cluster)[i])]] = c(mean(x), mean(y))
}
centers$`7`[1] = -14.5
centers$`7`[2] = 2
centers$`4`[1] = centers$`4`[1]-5
centers$`5`[1] = centers$`5`[1]+2
centers$`5`[2] = centers$`5`[2]-1

centers = data.frame(x = unlist(lapply(centers, function(x) x[1])),
                     y = unlist(lapply(centers, function(x) x[2])),
                     row.names = names(centers))

#Add to plot
for(i in 1:nrow(centers)){
  points(centers[i,1], 
         centers[i,2],
         pch = 21,
         cex = 3,
         bg = darken_color(RColorBrewer::brewer.pal(length(partition), 'Paired')[as.numeric(rownames(centers)[i])]))
}

#Occupancy across species
s = split(layout, layout$strain)
o = list()
for(i in 1:length(s)){
  x = split(s[[i]], s[[i]]$individual)
  o[[names(s)[i]]] = lapply(x, function(y) table(y$louvain_cluster)/nrow(y))
}

par(mfrow = c(2,4))
plot(layout[1:20000,1],
     layout[1:20000,2],
     pch = 20,
     col = alpha(cols, 0.5),
     cex = 0.25,
     bty = 'n',
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '')

for(i in 1:length(partition)){
  plot_results(lapply(o, function(x) unlist(lapply(x, function(y) y[i]))), 
               ylim = c(0, 1),
               ylab = '% time spent',
               xlab = '',
               col1 = darken_color(RColorBrewer::brewer.pal(length(partition), 'Paired')[i]),
               col2 = RColorBrewer::brewer.pal(length(partition), 'Paired')[i])
  
  k = kruskal.test(lapply(o, function(x) unlist(lapply(x, function(y) y[i]))))$p.value
  title(main = paste('p =', signif(k, 2)),
        font.main = 1,
        cex.main = 1.5)
}

#Vioplot
par(mfrow = c(2,4))
for(i in 1:length(o)){
  
  x1 = lapply(o, function(x) unlist(lapply(x, function(y) y[i])))
  x1 = lapply(x1, function(x) x[!is.na(x)])
  x = unlist(x1)
  g = unlist(lapply(strsplit(names(x), '_'), function(y) y[1]))
  
  p1 = signif(kruskal.test(x ~ as.factor(g))$p.value, 2)
  
  par(bty = 'n')
  stripchart(x ~ as.factor(g), 
             vertical = TRUE,
             pch = 20, 
             col = NULL,
             method = "jitter", 
             jitter = 0.2, 
             cex.lab = 1.5, 
             cex.axis = 1.5, 
             ylim = c(0,1),
             las = 2,
             xaxt = 'n',
             main = '',
             ylab = '% time spent')
  title(main = paste("p = ", p1),
        font.main = 1,
        cex.main = 1.5)
  
  vioplot.cols(x1,
               names = names(x1),
               col = rep(cols[i], length(x1)),
               border = rep(cols[i], length(x1)), 
               add = TRUE)
  
  for(i in 1:length(x1)){
    rect((i-1)+0.5, -2, (i-1)+1.5, 0, 
         col = dros_cols[names(dros_cols) == names(x1)[i]], 
         border = dros_cols[names(dros_cols) == names(x1)[i]])
  }
}

##Calculate transition probabilities
mc = markovchain::markovchainFit(layout$louvain_cluster)$estimate@transitionMatrix

#Chord diagram
cols = RColorBrewer::brewer.pal(length(partition), 'Paired')
cols[cols == "white"] = "grey90"
l_chord = circlize::chordDiagram(mc, 
                                 directional = TRUE,
                                 grid.col =  cols,
                                 self.link = 1,
                                 annotationTrack = c("name", "grid"))

par(mar = c(1,1,1,1))
plot_umap_markov(layout[1:10000,],
                 centers = centers,
                 l_chord,
                 plot_umap_points = FALSE,
                 cols = RColorBrewer::brewer.pal(nrow(centers), 'Paired')[as.numeric(rownames(centers))],
                 umap_cols = alpha(RColorBrewer::brewer.pal(nrow(centers), 'Paired')[layout[1:10000,]$louvain_cluster], 0.25))

#####################################################################################
#####Inter vs. intra-specific variation via Kruskal-Wallis test (as in Figure 3)#####
#####################################################################################
#Bin
layout = bin_umap(layout, n_bins = 100)$layout

#Function to calculate within species variance
within_species_variance_umap = function(layout,
                                        extract_condition = FALSE,
                                        condition){
  
  if(extract_condition == TRUE){
    #Extract species
    l = layout[layout$strain%in%condition,]
  }else{
    l = layout
  }
  
  #Split on trial
  trials = split(l, l$individual)
  
  #Extract all possible bins
  xy_new_u = unique(layout$xy_new)
  
  #Create dataframe for results
  n = names(trials)
  
  bin_variance = as.data.frame(matrix(nrow = length(xy_new_u),
                                      ncol = length(n)))
  colnames(bin_variance) = n
  rownames(bin_variance) = xy_new_u
  
  for(j in 1:length(trials)){
    #print(j)
    t = table(trials[[j]]$xy_new)
    
    bin_variance[,j] = t[match(rownames(bin_variance),
                               names(t))]
  }
  
  bin_variance[is.na(bin_variance)] = 0
  
  #Calculate variance
  v = apply(bin_variance, 1, function(x) sd(x))
  
  s = colSums(bin_variance)
  
  bin_variance_percent = bin_variance
  for(i in 1:ncol(bin_variance)){
    bin_variance_percent[,i] = bin_variance_percent[,i]/s[i]
  }
  
  v_p = apply(bin_variance_percent, 1, function(x) sd(x))
  
  z = list(bin_variance,
           v,
           bin_variance_percent,
           v_p)
  names(z) = c("bin_counts_raw",
               "bin_variance_raw",
               "bin_counts_percent",
               "bin_variance_percent")
  return(z)
}

#Run for all species
variance_maps = list()
for(h in 1:length(unique(layout$strain))){
  print(as.character(unique(layout$strain)[h]))
  variance_maps[[as.character(unique(layout$strain)[h])]] = within_species_variance_umap(layout,
                                                                                         extract_condition = TRUE,
                                                                                         unique(layout$strain)[h])$bin_counts_percent}

#Run kruskal.test
ps = c()
ks = c()
for(i in 1:nrow(variance_maps$la66)){
  print(i)
  ps = c(ps, kruskal.test(lapply(variance_maps[17:29], function(x) as.numeric(x[i,])))$p.value)
  ks = c(ks, kruskal.test(lapply(variance_maps[17:29], function(x) as.numeric(x[i,])))$statistic)
}

names(ps) = rownames(variance_maps$la66)

#Adjust
ps_a = ps*length(ps)

#Plot
ps_a = log10(ps)*-1
ps_a = round(ps_a, 2)
ps_a[ps_a<1.3] = 0

#Or statistic
ps_a = round(ks, 3)
names(ps_a) = rownames(variance_maps$la66)

#Colors
ints = seq(0, max(ps_a, na.rm = TRUE)+0.001, 0.001)
cols = colorRampPalette(c("white", "grey90", 'grey80', "red", 'darkred'))(length(ints))
names(cols) = ints
cols = cols[match(ps_a, names(cols))]
cols[is.na(cols)] = "grey90"

par(mar = c(2,2,2,2))
plot(unlist(lapply(strsplit(names(ps_a), "_"), function(v){v[1]})),
     unlist(lapply(strsplit(names(ps_a), "_"), function(v){v[2]})),
     pch = 20,
     cex = 0.75,
     #col = darken_color(cols),
     col = cols,
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     yaxt = 'n',
     ylab = "",
     xlab = "",
     bty = 'n')

#Legend
ints = seq(0, round(max(ks, na.rm = TRUE)), 1)
cols = colorRampPalette(c("white", "grey90", 'grey80', "red", 'darkred'))(length(ints))
names(cols) = ints

plot(seq(1,length(cols),1),
     rep(1, length(cols)),
     col = NULL,
     ylim = c(0,1),
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = 'Kruskal-Wallis Statistic',
     cex.lab = 1.5,
     bty = 'n')
for(i in 1:(length(cols)-1)){
  rect(i, 0, i+1, 1, col = cols[i], border = cols[i])
}
axis(1, 
     at = seq(0, 50, 10), 
     labels = seq(0, 50, 10),
     cex.axis = 1.5)

######################################################
#####Coverage of space by strains (as in Figure 3#####
######################################################
#Split on strain
strains = split(layout, as.character(layout$strain))

#Reorder by phylogeny
toMatch = c('simulans', 'sechellia', 'mauritiana', 'melanogaster', 'yakuba', 'santomea',
            'teissieri', 'erecta', 'persimilis', 'pseudoobscura', 'willistoni', 'virilis', 'arizonae',
            'la66', 'la69', 'z530', "z56", "z58", "zh16", "zh18", "zh20", "zh27",
            "zh29", "zh32", "zh33", "zh34", "zh42", "zh47", "zh58")
strains = strains[match(toMatch, names(strains))]

#Order colors
dros_cols = dros_cols[match(names(strains), names(dros_cols))]

#Get all unique xy coords in space
xy = unique(layout$xy_new)

props = list()
for(i in 1:length(xy)){
  print(i)
  props[[xy[i]]] = unlist(lapply(strains, function(x) nrow(x[x$xy_new == xy[i],])/nrow(x)))
}

mat = as.matrix(do.call(rbind, props))
mat = mat/rowSums(mat)
mat[is.na(mat)] = 0

d = dist(mat)
hcl = hclust(d)

par(mfrow = c(2,1), mar = c(1,1,1,1))
plot.phylo(as.phylo(hcl), direction = 'downwards', show.tip.label = FALSE)
barplot(t(mat[hcl$order,]), 
        col = dros_cols, 
        border = dros_cols,
        xaxt = 'n',
        yaxt = 'n')

#Overall coverage
cov = ncol(mat)-apply(mat, 1, function(x) sum(x==0))

#Histogram
par(mfrow = c(1,2))
cols = colorRamps::matlab.like(25)

h = hist(cov, 
         breaks = ncol(mat),
         xlab = 'n strains',
         cex.lab = 1.5,
         cex.axis = 1.5,
         las = 2,
         main = '',
         col = cols,
         border = darken_color(cols))
text(h$mids, 
     h$counts, 
     labels=round(h$counts, 3), 
     adj=c(0.5, -0.5))

#Coverage on space
cols = colorRamps::matlab.like(ncol(mat))
names(cols) = seq(1, ncol(mat), 1)
cols = cols[match(cov, names(cols))]

plot(unlist(lapply(strsplit(names(cov), "_"), function(v){v[1]})),
     unlist(lapply(strsplit(names(cov), "_"), function(v){v[2]})),
     pch = 21,
     cex = 1,
     bg = cols[match(cov, names(cols))],
     col = darken_color(cols[match(cov, names(cols))]),
     cex.axis = 1.5,
     cex.lab = 1.5,
     xaxt = 'n',
     yaxt = 'n',
     ylab = "",
     xlab = "",
     bty = 'n')

#########################################################
#####Behavior space autocorrelation (as in Figure 3)#####
#########################################################
#Calculate acfs
acfs = lapply(individuals, function(x) acf(x$coords, lag = 90, plot = FALSE)$acf[,,1])

#Calculate mean
m = rowMeans(do.call(cbind, acfs))

n = unlist(lapply(strsplit(names(acfs), "_"), function(v){v[1]}))
cols = dros_cols[match(n, names(dros_cols))]

#Plot
plot(m, 
     type = 'n',
     bty = 'n',
     ylab = "Correlation coefficient",
     xlab = 'Time (seconds)',
     cex.lab = 1.5,
     cex.axis = 1.5,
     xaxt = 'n',
     ylim = c(-0.5,1))
for(i in 1:length(acfs)){
  lines(acfs[[i]],
        col = alpha(cols[i], 0.5),
        lwd = 0.25)
}
lines(m,
      lwd = 3, 
      col = 'black')
abline(h = 0,
       lwd = 1.5,
       lty = 'dashed')
axis(1, 
     at = c(0, 30, 60, 90), 
     labels = c(0, 1, 2, 3),
     cex.axis = 1.5,
     cex.lab = 1.5)

######################################################
#####Calculate space persistence (as in Figure 3)#####
######################################################
#By louvain cluster
res = list()
for(h in 1:length(strains)){
  print(paste(h, 'out of', length(strains)))
  pos = unique(strains[[h]]$louvain_cluster)
  
  d = c()
  for(i in 1:length(pos)){
    x = which(strains[[h]]$louvain_cluster == pos[i])
    x = diff(x)
    d = c(d, x)
  }
  res[[names(strains)[h]]] = d
}

#Plot
hs = list()
for(i in 1:length(res)){
  
  z = res[[i]]
  z = z[z>1]
  z = z[z<90]
  hs[[names(strains)[i]]] = hist(z,
                                 xlim = c(0,90), 
                                 breaks = seq(1, 90, 1))$density
  hs[[names(strains)[i]]] = hs[[names(strains)[i]]]/sum(hs[[names(strains)[i]]])
}

hs = hs[-grep('^z', names(hs))]
hs = hs[-grep('la6', names(hs))]

plot(hs[[1]], 
     type = 'l', 
     col = dros_cols[match(names(hs)[1], names(dros_cols))],
     ylim = c(0, 0.12), 
     ylab = 'Porbability',
     xlab = 'Time (sec)', 
     cex.axis = 1.5, 
     cex.lab = 1.5, 
     bty = 'n',
     xaxt = 'n')
axis(1, c(0, 30, 60, 90), c(0, 1, 2, 3), cex.axis = 1.5)
for(i in 1:length(hs)){
  lines(hs[[i]], col = dros_cols[match(names(hs)[i], names(dros_cols))])
}

