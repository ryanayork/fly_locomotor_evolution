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

###################
#####Load data#####
###################
#Load layout
layout = readRDS("all_trials_umap_layout_annotated_30hz_size10.RDS")

#Load windows
win = readRDS('all_trials_velocity_windows_30hz_size10.RDS')

#Get behavior spaces by individual
individuals = split(layout, as.character(layout$individual))

#Get behavior spaces by strain
strains = split(layout, as.character(layout$strain))

##########################################################
#####Trial permutation power tests (as in Figure S2a)#####
##########################################################
#Calculate trial permutation correlations to full
strain_power_test = list()
for(a in 1:length(strains)){
  
  print(names(strains)[a])
  
  #Limits for pdfs
  xmin = min(layout$x)
  xmax = max(layout$x)
  ymin = min(layout$y)
  ymax = max(layout$y)
  
  #Get individual trials
  inds = split(strains[[a]], strains[[a]]$individual)
  
  #Get strain pdf
  s_pdf = unlist(as.data.frame(kde2d(strains[[a]]$x, strains[[a]]$y, h = 1, n = 100, lims = c(c(xmin, xmax), c(ymin, ymax)))$z))
  
  res = list()
  if(length(inds)>=10){
    for(i in 1:10){
      print(paste('n trials =', i))
      cors = c()
      for(j in 1:10){
        #Get random trials
        m = do.call(rbind, inds[sample(seq(1, length(inds), 1), i)])
        
        #Get trial pdf
        i_pdf = unlist(as.data.frame(kde2d(m$x, m$y, h = 1, n = 100, lims = c(c(xmin, xmax), c(ymin, ymax)))$z))
        
        #Correlate
        cors = c(cors, cor(s_pdf, i_pdf))
      }
      res[[as.character(i)]] = cors
    }
  }else{
    for(i in 1:length(inds)){
      print(paste('n trials =', i))
      cors = c()
      for(j in 1:length(inds)){
        #Get random trials
        m = do.call(rbind, inds[sample(seq(1, length(inds), 1), i)])
        
        #Get trial pdf
        i_pdf = unlist(as.data.frame(kde2d(m$x, m$y, h = 1, n = 100, lims = c(c(xmin, xmax), c(ymin, ymax)))$z))
        
        #Correlate
        cors = c(cors, cor(s_pdf, i_pdf))
      }
      res[[as.character(i)]] = cors
    }
  }
  
  plot_results(res, ylim = c(0,1), ylab = 'Correlation', xlab = 'n trials')
  
  strain_power_test[[names(strains)[a]]] = res
}

#Plot
par(mfrow = c(5,6), mar = c(2,2,2,2))
for(i in 1:length(strain_power_test)){
  plot_results(strain_power_test[[i]], 
               ylim = c(0,1), 
               ylab = 'Correlation', 
               xlab = 'n trials',
               col1 = darken_color(dros_cols[names(dros_cols)%in%names(strains[i])]),
               col2 = dros_cols[names(dros_cols)%in%names(strains[i])])
  title(main = names(strains)[i],
        cex.main = 1.5,
        font.main = 3,
        col.main = dros_cols[names(dros_cols)%in%names(strains[i])])}

#Plot all together
plot(unlist(lapply(strain_power_test[[1]], function(x) mean(x))),
     ylim = c(0,1),
     type = 'l',
     ylab = 'Correlation',
     xlab = 'n trials',
     col = alpha(dros_cols[names(dros_cols)%in%names(strains[1])], 0.5),
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5)
for(i in 2:length(strains)){
  lines(unlist(lapply(strain_power_test[[i]], function(x) mean(x))),
        lwd = 1.5,
        col = alpha(dros_cols[names(dros_cols)%in%names(strains[i])], 0.5))
}


###################################################################################################
#####Compare intra-trial, intra-species, and inter-species cosine distances (as in Figure S2b)#####
###################################################################################################
##Setup test sets
#Individuals
individuals2 = individuals[lapply(individuals, function(x) nrow(x))>5000]

#Strains (intra)
x = do.call(rbind, individuals2)
strains2 = split(x, x$strain)

##Run
#Individuals
ind_spec = list()
for(a in 1:length(individuals2)){
  
  print(paste(a, 'out of', length(individuals2)))
  
  #Limits for pdfs
  xmin = min(layout$x)
  xmax = max(layout$x)
  ymin = min(layout$y)
  ymax = max(layout$y)
  
  #Bootstrap
  res = list()
  for(i in 1:10){
    s = individuals2[[a]][sample(seq(1, nrow(individuals2[[a]]), 1), 2000),1:2]
    res[[i]] = unlist(as.data.frame(kde2d(s$x, s$y, h = 1, n = 100, lims = c(c(xmin, xmax), c(ymin, ymax)))$z))
  }
  res = do.call(cbind, res)
  
  #Calculate cosine distances in matrix
  d = lsa::cosine(res)
  d = unlist(as.data.frame(d))
  
  #Save
  ind_spec[[as.character(names(individuals2)[a])]] = d
}

#Strains intra
strains_spec = list()
for(a in 1:length(strains2)){
  
  print(paste(a, 'out of', length(strains2)))
  
  #Limits for pdfs
  xmin = min(layout$x)
  xmax = max(layout$x)
  ymin = min(layout$y)
  ymax = max(layout$y)
  
  #Bootstrap
  res = list()
  for(i in 1:100){
    s = strains2[[a]][sample(seq(1, nrow(strains2[[a]]), 1), 2000),1:2]
    res[[i]] = unlist(as.data.frame(kde2d(s$x, s$y, h = 1, n = 100, lims = c(c(xmin, xmax), c(ymin, ymax)))$z))
  }
  res = do.call(cbind, res)
  
  #Calculate cosine distances in matrix
  d = lsa::cosine(res)
  d = unlist(as.data.frame(d))
  
  #Save
  strains_spec[[as.character(names(strains2)[a])]] = d
}

#Strains inter
res = list()
for(a in 1:1000){
  
  if(a%%100 == TRUE){
    print(paste(a, 'out of', 1000))
  }
  
  #Limits for pdfs
  xmin = min(layout$x)
  xmax = max(layout$x)
  ymin = min(layout$y)
  ymax = max(layout$y)
  
  s = x[sample(seq(1, nrow(x), 1), 2000),1:2]
  res[[a]] = unlist(as.data.frame(kde2d(s$x, s$y, h = 1, n = 100, lims = c(c(xmin, xmax), c(ymin, ymax)))$z))
}

res = do.call(cbind, res)
strains_inter_spec = lsa::cosine(res)
strains_inter_spec = unlist(as.data.frame(strains_inter_spec))

##Compare
all = list(unlist(ind_spec), unlist(strains_spec), unlist(strains_inter_spec))
par(bty = 'n')
vioplot::vioplot(all[[1]][all[[1]]<1], 
                 all[[2]][all[[2]]<1],
                 all[[3]][all[[3]]<1][1:30000],
                 cex.axis = 1.5,
                 cex.lab = 1.5,
                 cex.sub = 1.5,
                 ylab = 'Cosine similarity',
                 las = 2,
                 ylim = c(0.7,1),
                 names = c('Individuals', 'Within species', 'Between species'),
                 col = 'grey80',
                 border = 'grey60')
title(main = paste('p =', signif(kruskal.test(all)$p.value, digits = 3)),
      cex.main = 1.5,
      font.main = 1.5)

#Post hoc
#Post hoc test
x = unlist(all)
g = c(rep('Individuals', length(all[[1]])),
      rep('Within species', length(all[[2]])),
      rep('Between species', length(all[[3]])))
d = dunn.test::dunn.test(x, g)

###############################################################################################################
#####Stationarity by comparing species pdfs made from the first half to the second half (as in Figure S2c)#####
###############################################################################################################
#Split individuals into first half and second half
first = list()
second = list()
for(i in 1:length(individuals)){
  
  n = nrow(individuals[[i]])
  n = round(n/2)
  
  first[[i]] = individuals[[i]][1:n,]
  second[[i]] = individuals[[i]][(n+1):nrow(individuals[[i]]),]
  
}

#Combine
first = do.call(rbind, first)
second = do.call(rbind, second)

#Split on species
first = split(first, first$strain)
second = split(second, second$strain)

#Caculate pdfs
first_pdfs = list()
second_pdfs = list()
for(i in 1:length(first)){
  print(paste(i, 'out of', length(first)))
  first_pdfs[[names(first)[i]]] = as.data.frame(kde2d(first[[i]]$x, first[[i]]$y, h = 1, n = 100)$z)
}
for(i in 1:length(second)){
  print(paste(i, 'out of', length(second)))
  second_pdfs[[names(second)[i]]] = as.data.frame(kde2d(second[[i]]$x, second[[i]]$y, h = 1, n = 100)$z)
}

#Correlate
cors = c()
for(i in 1:length(first)){
  cors = c(cors, cor(unlist(first_pdfs[[i]]), unlist(second_pdfs[[i]])))
}
names(cors) = names(first)
plot(cors, ylim = c(0,1))

#Just species
first_pdfs_s = first_pdfs[-grep('^z', names(first_pdfs))]
first_pdfs_s = first_pdfs_s[-grep('la6', names(first_pdfs_s))]

second_pdfs_s = second_pdfs[-grep('^z', names(second_pdfs))]
second_pdfs_s = second_pdfs_s[-grep('la6', names(second_pdfs_s))]

#Correlate
cors = c()
for(i in 1:length(first_pdfs_s)){
  cors = c(cors, cor(unlist(first_pdfs_s[[i]]), unlist(second_pdfs_s[[i]])))
}
names(cors) = names(first_pdfs_s)
plot(cors, 
     ylim = c(0,1),
     ylab = 'Correlation',
     xlab = '',
     cex.axis = 1.5, 
     las = 2,
     cex.lab = 1.5, bty = 'n', pch = 20, cex = 2, 
     col = dros_cols[match(names(cors), names(dros_cols))], xaxt = 'n')
axis(1, 1:length(cors), names(cors), las = 2, cex.lab = 1.5, cex.axis = 1.5)

#Plot
par(mfrow = c(2, 13), mar = c(1,1,1,1))
for(i in 1:length(first_pdfs_s)){
  image(as.matrix(first_pdfs_s[[i]]), col = colorRamps::matlab.like(100), xaxt = 'n', yaxt = 'n', bty = 'n')
  title(main = names(first_pdfs_s)[i], cex.main = 1.5, font.main = 1)
}
for(i in 1:length(second_pdfs_s)){
  image(as.matrix(second_pdfs_s[[i]]), col = colorRamps::matlab.like(100), xaxt = 'n', yaxt = 'n', bty = 'n')
}

