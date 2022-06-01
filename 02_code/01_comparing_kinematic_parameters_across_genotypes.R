############################
#####Initiate workspace#####
############################
rm(list=ls());
options(stringsAsFactors=F);
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

####################################################################################
#####Compare velocity measures as a function of strain (as in Figures 1 and S1)#####
####################################################################################
#Load data
vel = readRDS('dros_locomotor_evolution_all_velocity_smoothed_xy_092520.RDS')

#Overall autocorrelations
acf(as.numeric(na.omit(unlist(lapply(vel, function(x) x$vr_cm_smooth)))), lag = 300)

#Get species-wise vr
vr = unlist(lapply(vel, function(x) x$vr_cm_smooth))
s = unlist(lapply(vel, function(x) x$strain))

vr = split(vr, s)

#Calculate species-wise acf
acfs = lapply(vr, function(x) acf(as.numeric(na.omit(x)), lag = 100, plot = FALSE)$acf[,,1])

#Plot
par(mfrow = c(3,8), mar = c(2,2,2,2))
for(i in 1:length(acfs)){
  plot(acfs[[i]],
       ylim = c(-1,1),
       type = 'h',
       xlab = 'Time',
       ylab = 'Correlation',
       cex.axis = 1.5,
       cex.lab = 1.5,
       bty = 'n',
       col = 'gray70',
       lwd = 3)
  abline(h = 0, lwd = 1.5)
  abline(v = which.min(acfs[[i]]>0),
         col = 'red',
         lty = 'dashed', 
         lwd = 1.5)
  title(main = paste(names(acfs)[i], '=', which.min(acfs[[i]]>0)),
        cex.main = 1.5,
        font.main = 1)
}

#Velocity by strain
vt = unlist(lapply(vel, function(x) x$vt_cm))
s = unlist(lapply(vel, function(x) x$strain))

vt = split(vt, s)

par(mfrow = c(5,6))
for(i in 1:length(vt)){
  hist(vt[[i]], 
       xlim = c(0,3.5),
       main = '',
       xlab = 'Velocity (cm/sec)',
       cex.axis = 1.5,
       cex.lab = 1.5,
       breaks = 50,
       col = dros_cols[names(dros_cols)==names(vt)[i]],
       border = dros_cols[names(dros_cols)==names(vt)[i]])
  title(main = names(vt)[i],
        cex.main = 1.5,
        font.main = 3,
        col.main = dros_cols[names(dros_cols)==names(vt)[i]])
}

#Overall velocities
vt = unlist(lapply(vel, function(x) x$vt_cm))
vr = unlist(lapply(vel, function(x) x$vr_cm))

par(mfrow = c(1,3))
hist(vt, 
     xlim = c(0,3), 
     breaks = 100,
     cex.lab = 1.5,
     cex.axis = 1.5,
     col = 'grey70',
     xlab = 'Translational velocity (cm/sec)',
     main = '',
     border = 'grey70')

hist(vr,
     breaks = 200,
     cex.lab = 1.5,
     cex.axis = 1.5,
     col = 'grey70',
     xlab = 'Angular velocity (deg/sec)',
     main = '',
     border = 'grey70')

#Distance traveled
d = unlist(lapply(vel, function(x) sum(pracma::hypot(diff(x$x), diff(x$y)))))
hist(d,
     breaks = 50,
     cex.lab = 1.5,
     cex.axis = 1.5,
     col = 'grey70',
     xlab = 'Distance covered',
     main = '',
     border = 'grey70')

#################################################################################
#####Analyzing distribution of time spent moving by strain (as in Figure S1)#####
#################################################################################
#Load velocity data
vel = readRDS('~/Desktop/fly_locomotor_evolution_ms/00_data/00_velocity_files/dros_locomotor_evolution_all_velocity_smoothed_xy_092520.RDS')

#Loop through trials and calculate %of time moving (>0.05 cm/sec)
move = lapply(vel, function(x) sum(na.omit(x$vt_cm)>0.05)/nrow(x))
move = unlist(move)*100

#Overall histogram
h <- hist(move,
          breaks = 10, 
          plot = FALSE)
h$counts=h$counts/sum(h$counts)
plot(h,
     cex.lab = 1.5,
     cex.axis = 1.5,
     col = 'grey70',
     xlab = '% of time spent moving',
     ylab = 'Probability',
     main = '',
     border = 'grey70')

#Split on strain
move = split(move, unlist(lapply(strsplit(names(move), "_"), function(v){v[1]})))
move = move[match(names(dros_cols), names(move))]

#Plot
plot(1:length(move),
     rep(1, length(move)),
     ylim = c(0, 100),
     xlim = c(1, length(move)),
     type = 'n',
     ylab = '% of time spent moving',
     xlab = '',
     xaxt = 'n',
     bty = 'n', 
     las = 2,
     cex.axis = 1.5,
     cex.lab = 1.5)
for(j in 1:length(move)){
  s = boxplot.stats(unlist(move[[j]]))
  points(j+0.25, s$stats[3], pch = 20, cex = 3, col = darken_color(dros_cols[j]))
  segments(j+0.25, s$stats[2], j+0.25, s$stats[4], lwd = 1.5, col = darken_color(dros_cols[j]))
}
axis(1, at = 1.25:(length(move)+0.25), labels = names(move),
     cex.axis = 1.5,
     cex.lab = 1.5,
     las = 2)



