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

##################################
#####Plot phylogeny with time#####
##################################
#Load tree
tree = read.tree('fossilized_filtered_tree.tre')

#Plot
plotTree(tree, 
         direction="leftwards", 
         xlim=c(50,-15), 
         mar = c(4,2,2,4), 
         lwd = 1)
axis(1, at=seq(0, 45, by=5))

#########################
#####Load trait data#####
#########################
#Load layout
layout = readRDS('all_trials_umap_layout_annotated_30hz_size10_with_louvain_clusters_and_position.RDS')

#Split into strains
strains = split(layout, layout$strain)

#Split on individual
individuals = split(layout, as.character(layout$individual))

#############################
#####Calculate frequency#####
#############################
#Get all unique xy coords in space
xy = unique(layout$xy_new)

#Calculate space coverage by individual
percs = list()

for(i in 1:length(individuals)){
  x = unique(individuals[[i]]$xy_new)
  percs[[names(individuals)[i]]] = length(x%in%xy)/length(xy)
}
names(percs) = unlist(lapply(strsplit(names(percs), "_"), function(v){v[1]}))

#Calculate pdfs
pdfs = list()
xmin = min(layout$x)
xmax = max(layout$x)
ymin = min(layout$y)
ymax = max(layout$y)

for(i in 1:length(individuals)){
  print(paste(i, 'out of', length(individuals)))
  pdfs[[names(individuals)[i]]] = kde2d(individuals[[i]]$x,
                                        individuals[[i]]$y,
                                        h = 1,
                                        n = 100,
                                        #n = 32,
                                        lims = c(c(xmin, xmax),
                                                 c(ymin, ymax)))}

#Unlist and combine into matrix
ind_pdfs = do.call(cbind, lapply(pdfs, function(x) unlist(as.data.frame(x$z))))
ind_pdfs = apply(ind_pdfs, 2, function(x) x/max(x))

#Calculate mean and se
s = unique(unlist(lapply(strsplit(colnames(ind_pdfs), '_'), function(y) y[1])))

pdf_means = list()
pdf_ses = list()
for(i in 1:length(s)){
  
  x = ind_pdfs[,grep(s[i], colnames(ind_pdfs))]
  
  pdf_means[[s[i]]] = rowMeans(x)
  pdf_ses[[s[i]]] =  apply(x, 1, function(y) plotrix::std.error(y))
}

pdf_means = do.call(rbind, pdf_means)
pdf_ses = do.call(rbind, pdf_ses)

pdf_means = pdf_means[match(tree$tip.label, rownames(pdf_means)),]
pdf_ses = pdf_ses[match(tree$tip.label, rownames(pdf_ses)),]

###############################
#####Calculate transitions#####
###############################
#Calculate markov models and extract transition probabilities
partition = unique(layout$louvain_cluster)
chords = list()
mcs = list()
for(i in 1:length(strains)){
  
  ##Calculate transition probabilities
  mc = markovchain::markovchainFit(strains[[i]]$louvain_cluster)
  
  m = mc$estimate@transitionMatrix
  
  #Chord diagram
  cols = RColorBrewer::brewer.pal(length(partition), 'Paired')
  cols[cols == "white"] = "grey90"
  chords[[names(strains)[i]]] = circlize::chordDiagram(m, 
                                                       directional = TRUE,
                                                       direction.type = "arrows",
                                                       grid.col =  cols,
                                                       self.link = 1,
                                                       annotationTrack = c("name", "grid"))
  
  mcs[[names(strains)[i]]] = mc
}

#Combine into matrix
trait = data.frame(matrix(nrow = length(tree$tip.label), ncol = 0))
for(i in 1:length(unlist(as.data.frame(mcs$arizonae$estimate@transitionMatrix)))){
  x = unlist(lapply(mcs, function(y) unlist(as.data.frame(y$estimate@transitionMatrix))[i]))
  x = x[match(tree$tip.label, 
              unlist(lapply(strsplit(names(x), '\\.'), function(y) y[1])))]
  trait = cbind(trait, x)
}
colnames(trait) = names(unlist(as.data.frame(mcs$arizonae$estimate@transitionMatrix)))
rownames(trait) = tree$tip.label

#Variation in markov models across strains factoring in individuals
strain_mcs = list()
for(i in 1:length(strains)){
  
  print(i)
  inds = split(strains[[i]], strains[[i]]$individual)
  
  tmp = list()
  for(j in 1:length(inds)){
    
    z = markovchain::markovchainFit(inds[[j]]$louvain_cluster)$estimate@transitionMatrix
    
    z = z[,match(seq(1, 7, 1), colnames(z))]
    z = z[match(seq(1, 7, 1), rownames(z)),]
    z[is.na(z)] = 0
    
    rownames(z) = seq(1,7,1)
    colnames(z) = seq(1,7,1)
    
    tmp[[names(inds)[j]]] = unlist(as.data.frame(z))
    
  }
  strain_mcs[[names(strains)[i]]] = do.call(cbind, tmp)
}

#Calculate mean and se
means = do.call(rbind, lapply(strain_mcs, function(x) rowMeans(x)))
ses = do.call(rbind, lapply(strain_mcs, function(x) apply(x, 1, function(y) plotrix::std.error(y))))

means = means[match(tree$tip.label, rownames(means)),]
ses = ses[match(tree$tip.label, rownames(ses)),]

#################################################
#####Calculate variation in louvain clusters#####
#################################################
#Calculate cluster occupancy by strain
clusters = do.call(rbind, lapply(strains, function(x) table(x$louvain_cluster)/nrow(x)))

#Match to tree
clusters = clusters[match(tree$tip.label, rownames(clusters)),]

#############################
#####Calculate structure#####
#############################
#Calculate binary occupancy per strain
xy = unique(layout$xy_new)
binary = list()
for(i in 1:length(strains)){
  print(i)
  
  x = xy%in%unique(strains[[i]]$xy_new)
  
  x[x == TRUE] = 1
  x[!x==1] = 0
  names(x) = xy
  binary[[names(strains)[i]]] = x
}

#Combine
binary = do.call(rbind, binary)

#Match to tree
binary = binary[match(tree$tip.label, rownames(binary)),]

##################################################################################################
#####Generating and plotting binary, frequency, and transitions morphospaces (as in Figure 5)#####
##################################################################################################
#Binary
pca = prcomp(cbind(binary))

plot(pca$x[,1:2], 
     pch = 20,
     col = dros_cols[match(rownames(pca$x), names(dros_cols))],
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlab = paste('PC1', '(30.85% variation explained)'),
     ylab = paste('PC2', '(12.47% variation explained)'))
text(pca$x[,1:2], 
     rownames(pca$x),
     col = dros_cols[match(rownames(pca$x), names(dros_cols))])

#Occupancy
pca = prcomp(cbind(pdf_means))

plot(pca$x[,1:2], 
     pch = 20,
     col = dros_cols[match(rownames(pca$x), names(dros_cols))],
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlab = paste('PC1', '(51.73% variation explained)'),
     ylab = paste('PC2', '(20.23% variation explained)'))
text(pca$x[,1:2], 
     rownames(pca$x),
     col = dros_cols[match(rownames(pca$x), names(dros_cols))])

#Density maps of first two PCs
par(mfrow = c(1,2))
x = abs(pca$rotation[,1])
x[x<0.00001] = 0
image(matrix(x, nrow = 100),
      col = colorRampPalette(c('white', 'grey70', 'black'))(100),
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n')

x = abs(pca$rotation[,2])
x[x<0.00001] = 0
image(matrix(x, nrow = 100),
      col = colorRampPalette(c('white', 'grey70', 'black'))(100),
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n')

#Markov
pca = prcomp(cbind(means))

plot(pca$x[,1:2], 
     pch = 20,
     col = dros_cols[match(rownames(pca$x), names(dros_cols))],
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlab = paste('PC1', '(66.97% variation explained)'),
     ylab = paste('PC2', '(13.75% variation explained)'))
text(pca$x[,1:2], 
     rownames(pca$x),
     col = dros_cols[match(rownames(pca$x), names(dros_cols))])

################################################################
#####Plot eigenvectors for all morphospaces (as in Figure 5#####
################################################################
##Transitions
#Run pca
pca = prcomp(cbind(means))

#Generate plot parameters
res = list()
for(h in 1:3){
  x = pca$rotation[,h]
  x = round(x, 2)
  
  #Calculate centers
  partition = unique(layout$louvain_cluster)
  
  #Calculate centers
  centers = list()
  for(i in 1:length(partition)){
    x1 = layout$x[layout$louvain_cluster == unique(layout$louvain_cluster)[i]]
    y = layout$y[layout$louvain_cluster == unique(layout$louvain_cluster)[i]]
    centers[[as.character(unique(layout$louvain_cluster)[i])]] = c(mean(x1), mean(y))
  }
  centers$`7`[1] = -14.5
  centers$`7`[2] = 2
  centers$`4`[1] = centers$`4`[1]-5
  centers$`5`[1] = centers$`5`[1]+2
  centers$`5`[2] = centers$`5`[2]-1
  
  centers = data.frame(x = unlist(lapply(centers, function(x) x[1])),
                       y = unlist(lapply(centers, function(x) x[2])),
                       row.names = names(centers))
  
  #Variation in transitions
  mc_var = as.data.frame(matrix(x, nrow = 7))
  
  ##Calculate transition probabilities
  mc = markovchain::markovchainFit(layout$louvain_cluster)$estimate@transitionMatrix
  
  #Chord diagram
  cols = RColorBrewer::brewer.pal(length(partition), 'Paired')
  cols[cols == "white"] = "grey90"
  chord = circlize::chordDiagram(mc, 
                                 directional = TRUE,
                                 grid.col =  cols,
                                 self.link = 1,
                                 annotationTrack = c("name", "grid"))
  
  #Get colors
  cols1 = colorRampPalette(c("midnightblue", "grey90"))(length(seq(min(x), 0, 0.01)))
  names(cols1) = seq(min(x), 0, 0.01)
  cols2 = colorRampPalette(c("grey90", 'darkred'))(length(seq(0.01, max(x), 0.01)))
  names(cols2) = seq(0.01, max(x), 0.01)
  
  cols = c(cols1, cols2)
  cols = cols[match(x, names(cols))]
  names(cols) = paste(chord$rn, chord$cn, sep = '')
  
  #Plotting in umap
  #Reformat chord diagram coords
  chord$fromx = centers[match(chord$rn, rownames(centers)), 1]
  chord$fromy = centers[match(chord$rn, rownames(centers)), 2]
  chord$tox = centers[match(chord$cn, rownames(centers)), 1]
  chord$toy = centers[match(chord$cn, rownames(centers)), 2]
  
  #Split into self and not
  chord_s = chord[chord$rn == chord$cn,]
  chord = chord[!chord$rn == chord$cn,]
  
  #Make colors not alpha
  chord_s$col = cols[match(paste(chord_s$rn, chord_s$cn, sep = ''),
                           names(cols))]
  chord$col = cols[match(paste(chord$rn, chord$cn, sep = ''),
                         names(cols))]
  
  #Remove zeroes
  chord = chord[!chord$value1==0,]
  
  l = list(chord, chord_s, centers)
  names(l) = c('chord', 'chord_s', 'centers')
  res[[as.character(h)]] = l
}

#Plot
par(mfrow = c(1,3))
for(i in 1:3){
  
  chord = res[[i]]$chord
  chord_s = res[[i]]$chord_s
  centers = res[[i]]$centers
  
  plot(layout$x[1:1000],
       layout$y[1:1000],
       pch = 20,
       #col = "grey90",
       col = NULL,
       bty = 'n',
       xaxt = 'n',
       yaxt = 'n',
       xlab = "",
       ylab = "")
  
  for(i in 1:nrow(chord)){
    diagram::curvedarrow(from = c(chord$fromx[i], chord$fromy[i]),
                         to = c(chord$tox[i], chord$toy[i]),
                         lcol = chord$col[i],
                         lwd = chord$value1[i]*20,
                         arr.pos = 0.9,
                         curve = 0.2,
                         arr.type = "triangle",
                         arr.length = 0,
                         arr.width = 0)}
  for(i in 1:nrow(chord)){
    diagram::selfarrow(c(chord_s$fromx[i], chord_s$fromy[i]),
                       lcol = chord_s$col[i],
                       lwd = chord_s$value1[i]*10,
                       curve = 1,
                       arr.length = 0,
                       arr.width = 0)}
  
  points(centers,
         pch = 21,
         bg = RColorBrewer::brewer.pal(nrow(centers), 'Paired')[as.numeric(rownames(centers))],
         col = 'black',
         cex = 3)
}

##Binary
#Run pca
pca = prcomp(cbind(binary))

#Plot first n PCs
#Get all possible bins (64x64 grid)
all_bins = expand.grid(1:64, 1:64)
all_bins = paste(all_bins[,1], all_bins[,2], sep = '_')

#Plot n PCs
par(mfrow = c(3,3))
for(i in 1:3){
  x = pca$rotation[,i]
  x = round(x, 2)
  
  cols1 = colorRampPalette(c("midnightblue", "grey90"))(length(seq(min(x), 0, 0.01)))
  names(cols1) = seq(min(x), 0, 0.01)
  cols2 = colorRampPalette(c("grey90", 'darkred'))(length(seq(0.01, max(x), 0.01)))
  names(cols2) = seq(0.01, max(x), 0.01)
  
  cols = c(cols1, cols2)
  cols = cols[match(x, names(cols))]
  
  plot(unlist(lapply(strsplit(rownames(pca$rotation), "_"), function(v){v[1]})),
       unlist(lapply(strsplit(rownames(pca$rotation), "_"), function(v){v[2]})),
       col = cols,
       pch = 20,
       bty = 'n',
       xaxt = 'n',
       yaxt = 'n', 
       ylab = '', 
       xlab = '')
}

##Frequency
#Run pca
pca = prcomp(cbind(pdf_means))

#Plot n pcs
#par(mfrow = c(1,3))
for(i in 1:3){
  x = pca$rotation[,i]
  
  cols = make_centered_colors(minimum = min(x), maximum = max(x), 'darkblue', 'white', 'darkred', seq_by = 0.001)
  image(matrix(x, nrow = 100), col = cols$rampcols, breaks = cols$rampbreaks, bty = 'n', xaxt = 'n', yaxt = 'n')
}

#Add transitions
for(i in 1:3){
  
  chord = res[[i]]$chord
  chord_s = res[[i]]$chord_s
  centers = res[[i]]$centers
  
  plot(layout$x[1:1000],
       layout$y[1:1000],
       pch = 20,
       #col = "grey90",
       col = NULL,
       bty = 'n',
       xaxt = 'n',
       yaxt = 'n',
       xlab = "",
       ylab = "")
  
  for(i in 1:nrow(chord)){
    diagram::curvedarrow(from = c(chord$fromx[i], chord$fromy[i]),
                         to = c(chord$tox[i], chord$toy[i]),
                         lcol = chord$col[i],
                         lwd = chord$value1[i]*20,
                         arr.pos = 0.9,
                         curve = 0.2,
                         arr.type = "triangle",
                         arr.length = 0,
                         arr.width = 0)}
  for(i in 1:nrow(chord)){
    diagram::selfarrow(c(chord_s$fromx[i], chord_s$fromy[i]),
                       lcol = chord_s$col[i],
                       lwd = chord_s$value1[i]*10,
                       curve = 1,
                       arr.length = 0,
                       arr.width = 0)}
  
  points(centers,
         pch = 21,
         bg = RColorBrewer::brewer.pal(nrow(centers), 'Paired')[as.numeric(rownames(centers))],
         col = 'black',
         cex = 3)
}

#########################################################################################################################
#####Phylogenetically independent contrasts (to get rate): Space occupancy by bin (rate landscape of behavior space)#####
#########################################################################################################################
#Calculate pics
pics = apply(pdf_means, 2, function(x) median(abs(pic(x, tree))))

#Filter
#pics[pics<quantile(pics, probs = 0.4)] = 0

#Plot
#image(matrix(pics, ncol = 32, nrow = 32), col = c('white', hcl.colors(12, "YlOrRd", rev = TRUE)))

#3d scatter
x = expand.grid(seq(1,100,1), seq(1,100,1))[,1]
y = expand.grid(seq(1,100,1), seq(1,100,1))[,2]
z = pics
#x = x[!z == 0]
#y = y[!z == 0]
#z = z[!z == 0]

library(plot3D)

mod <- mgcv::gam(z ~ s(x, y, k=150))
m2 <- matrix(fitted(mod), ncol = 100)

persp3D(z=m2,
        col = colorRampPalette(c('white', 'grey90', 'grey80',
                                 'lightgoldenrod1', "gold","tomato",
                                 "darkred", "#3D0404", 'cadetblue1'))(100),
        box = FALSE,
        axis = FALSE, 
        grid = FALSE,
        expand = 0.5, 
        contour =  list(side = c("z")))

#2d
image(matrix(unlist(pics), nrow = 100),
      col = colorRampPalette(c('white', 'grey90', 'grey80',
                               "red", "darkred"))(200),
      bty = 'n',
      xaxt = 'n',
      yaxt = 'n')

#Phylogenetic models
delta_lik = function(x) x-x[which(x==max(x))]

mods = list()
best = c()
for(i in 1:ncol(pdf_means)){
  
  print(i)
  bm = evo.model(tree = mytimetree,
                 Y = pdf_means[,i],
                 model = 'BM') 
  ou = evo.model(tree = mytimetree,
                 Y = pdf_means[,i],
                 model = 'OU')
  eb = evo.model(tree = mytimetree,
                 Y = pdf_means[,i],
                 model = 'EB')
  
  liks = unlist(lapply(list(mod_bm,
                            mod_ou,
                            mod_eb),
                       function(x) x$logL))
  names(liks) = c('bm', 'ou', 'eb')
  liks = delta_lik(liks)
  
  best = c(best, names(liks)[which.max(liks)])
  mods[[i]] = liks
}

#########################################################################################################
#####Phylogenetically independent contrasts (to get rate): Markov (rate landscape of behavior space)#####
#########################################################################################################
#Calculate pics
pics = apply(means, 2, function(x) median(abs(pic(x, tree))))

#Calculate centers
partition = unique(layout$louvain_cluster)

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

#Variation in transitions
mc_var = as.data.frame(matrix(pics, nrow = 7))
#chord = chords$arizonae
#chord$value1 = unlist(mc_var)
#chord$value2 = unlist(mc_var)

##Calculate transition probabilities
mc = markovchain::markovchainFit(layout$louvain_cluster)$estimate@transitionMatrix

#Chord diagram
cols = RColorBrewer::brewer.pal(length(partition), 'Paired')
cols[cols == "white"] = "grey90"
chord = circlize::chordDiagram(mc, 
                               directional = TRUE,
                               grid.col =  cols,
                               self.link = 1,
                               annotationTrack = c("name", "grid"))

#Get colors
cols = colorRampPalette(c('grey90', 'grey80',
                          "red", "darkred"))(length(seq(0, max(round(pics, 2)), 0.01)))

names(cols) = seq(0, max(round(pics, 2)), 0.01)
cols = cols[match(round(pics, 2), names(cols))]
names(cols) = paste(chord$rn, chord$cn, sep = '')

#Plotting in umap
#Reformat chord diagram coords
chord$fromx = centers[match(chord$rn, rownames(centers)), 1]
chord$fromy = centers[match(chord$rn, rownames(centers)), 2]
chord$tox = centers[match(chord$cn, rownames(centers)), 1]
chord$toy = centers[match(chord$cn, rownames(centers)), 2]

#Split into self and not
chord_s = chord[chord$rn == chord$cn,]
chord = chord[!chord$rn == chord$cn,]

#Make colors not alpha
chord_s$col = cols[match(paste(chord_s$rn, chord_s$cn, sep = ''),
                         names(cols))]
chord$col = cols[match(paste(chord$rn, chord$cn, sep = ''),
                       names(cols))]

#Remove zeroes
chord = chord[!chord$value1==0,]

plot(layout$x,
     layout$y,
     pch = 20,
     #col = "grey90",
     col = NULL,
     bty = 'n',
     xaxt = 'n',
     yaxt = 'n',
     xlab = "",
     ylab = "")


for(i in 1:nrow(chord)){
  diagram::curvedarrow(from = c(chord$fromx[i], chord$fromy[i]),
                       to = c(chord$tox[i], chord$toy[i]),
                       lcol = chord$col[i],
                       lwd = chord$value1[i]*20,
                       arr.pos = 0.9,
                       curve = 0.2,
                       arr.type = "triangle",
                       arr.length = 0,
                       arr.width = 0)}

for(i in 1:nrow(chord)){
  diagram::selfarrow(c(chord_s$fromx[i], chord_s$fromy[i]),
                     lcol = chord_s$col[i],
                     lwd = chord_s$value1[i]*10,
                     curve = 1,
                     arr.length = 0,
                     arr.width = 0)}

points(centers,
       pch = 21,
       bg = RColorBrewer::brewer.pal(nrow(centers), 'Paired')[as.numeric(rownames(centers))],
       col = 'black',
       cex = 3)

##################################################################################################
#####Calculate phylogenetic signal: Space occupancy by bin (binary by strain; i.e. structure)#####
##################################################################################################
#Calculate phylogenetic signal
binary_k = list()
for(i in 1:ncol(binary)){
  print(i)
  #out[[i]] = phylosig(tree, pdf_means[,i], se = pdf_ses[,i], nsim = 10000, test = TRUE)
  
  binary_k[[i]] = phylosig(tree, binary[,i])
  print(binary_k[[i]])
  #print(out[[i]]$P)
}

#Alternative: %space covered per species
xy = unique(layout$xy_new)
binary = list()
for(i in 1:length(strains)){
  print(i)
  
  x = xy%in%unique(strains[[i]]$xy_new)
  x = sum(x == TRUE)/length(x)
  binary[[names(strains)[i]]] = x
  
}

binary = do.call(rbind, binary)
binary = binary[match(tree$tip.label, rownames(binary)),]

###############################################################
#####Calculate phylogenetic signal: Space occupancy by bin#####
###############################################################
#Calculate phylogenetic signal factoring in standard error
out = list()
for(i in 1:ncol(pdf_means)){
  print(i)
  #out[[i]] = phylosig(tree, pdf_means[,i], se = pdf_ses[,i], nsim = 10000, test = TRUE)
  
  if(sum(pdf_means[,i]<quantile(unlist(pdf_means), prob = 0.3))>10){
    out[[i]] = list(K = NA)
  }else{
    out[[i]] = phylosig(tree, pdf_means[,i], se = pdf_ses[,i])
    print(out[[i]]$K)
  }
  #print(out[[i]]$P)
}

#Distribution of K
#occupancy_k = matrix(unlist(lapply(out, function(x) x$K)), ncol = 64, nrow = 64)
occupancy_k = matrix(unlist(lapply(out, function(x) x$K)), ncol = 100, nrow = 100)
hist(occupancy_k)

#Plot
image2D(occupancy_k, col = colorRampPalette(c('white', 'grey90', 'grey80',
                                              'lightgoldenrod1', "gold","tomato",
                                              "red", "darkred", "#3D0404"))(100))

########################################################################
#####Calculate phylogenetic signal: Markov transition probabilities#####
########################################################################
#Get n time points per strain
s = lapply(strains, function(x) x[sample(1:nrow(x), 20000, replace = FALSE),])
s = do.call(rbind, s)

#Loop through and test
ks = list()
for(i in 1:10){
  
  newlayout = bin_umap(s, n_bins = i)
  
  #Get all unique xy coords in space
  xy = unique(newlayout$xy_new)
  
  #Calculate space coverage by individual
  percs = list()
  
  for(i in 1:length(individuals)){
    x = unique(individuals[[i]]$xy_new)
    percs[[names(individuals)[i]]] = length(x%in%xy)/length(xy)
  }
  names(percs) = unlist(lapply(strsplit(names(percs), "_"), function(v){v[1]}))
  
  #Calculate pdfs
  pdfs = list()
  xmin = min(layout$x)
  xmax = max(layout$x)
  ymin = min(layout$y)
  ymax = max(layout$y)
  
  for(i in 1:length(individuals)){
    print(paste(i, 'out of', length(individuals)))
    pdfs[[names(individuals)[i]]] = kde2d(individuals[[i]]$x,
                                          individuals[[i]]$y,
                                          h = 1,
                                          n = 100,
                                          #n = 32,
                                          lims = c(c(xmin, xmax),
                                                   c(ymin, ymax)))}
  
  #Unlist and combine into matrix
  ind_pdfs = do.call(cbind, lapply(pdfs, function(x) unlist(as.data.frame(x$z))))
  ind_pdfs = apply(ind_pdfs, 2, function(x) x/max(x))
  
  #Calculate mean and se
  s = unique(unlist(lapply(strsplit(colnames(ind_pdfs), '_'), function(y) y[1])))
  
  pdf_means = list()
  pdf_ses = list()
  for(i in 1:length(s)){
    
    x = ind_pdfs[,grep(s[i], colnames(ind_pdfs))]
    
    pdf_means[[s[i]]] = rowMeans(x)
    pdf_ses[[s[i]]] =  apply(x, 1, function(y) plotrix::std.error(y))
  }
  
  pdf_means = do.call(rbind, pdf_means)
  pdf_ses = do.call(rbind, pdf_ses)
  
  pdf_means = pdf_means[match(tree$tip.label, rownames(pdf_means)),]
  pdf_ses = pdf_ses[match(tree$tip.label, rownames(pdf_ses)),]
  
}

#Phylogenetic signal
out = apply(trait, 2, function(x) multiPhylosignal(as.data.frame(x), tree, reps = 1000))
ps = unlist(lapply(out, function(x) x$PIC.variance.P))
k = unlist(lapply(out, function(x) x$K))

#Phylogenetic signal factoring in standard error
out = list()
for(i in 1:ncol(means)){
  print(i)
  #out[[colnames(means)[i]]] = phylosig(tree, means[,i], se = ses[,i], nsim = 10000, test = TRUE)
  out[[colnames(means)[i]]] = phylosig(tree, means[,i], se = ses[,i])
  print(out[[i]]$P)
}

#Distribution of K
markov_k = unlist(lapply(out, function(x) x$K))
hist(markov_k)

############################################################################################
#####Calculate phylogenetic signal: frequency as a function of n bins in behavior space#####
############################################################################################
#Get n time points per strain
newstrains = lapply(strains, function(x) x[sample(1:nrow(x), 20000, replace = FALSE),])
newstrains = do.call(rbind, newstrains)

ks = list()
for(h in 2:100){
  
  print(paste('bin size', h))
  
  #Calculate pdfs
  pdfs = list()
  xmin = min(layout$x)
  xmax = max(layout$x)
  ymin = min(layout$y)
  ymax = max(layout$y)
  
  print('getting pdfs')
  for(i in 1:length(individuals)){
    #print(paste(i, 'out of', length(newinds)))
    pdfs[[names(newinds)[i]]] = kde2d(individuals[[i]]$x,
                                      individuals[[i]]$y,
                                      h = 1,
                                      n = h,
                                      lims = c(c(xmin, xmax),
                                               c(ymin, ymax)))}
  
  #Unlist and combine into matrix
  ind_pdfs = do.call(cbind, lapply(pdfs, function(x) unlist(as.data.frame(x$z))))
  ind_pdfs = apply(ind_pdfs, 2, function(x) x/max(x))
  
  #Calculate mean and se
  s = unique(unlist(lapply(strsplit(colnames(ind_pdfs), '_'), function(y) y[1])))
  
  new_pdf_means = list()
  new_pdf_ses = list()
  for(i in 1:length(s)){
    
    x = ind_pdfs[,grep(s[i], colnames(ind_pdfs))]
    
    new_pdf_means[[s[i]]] = rowMeans(x)
    new_pdf_ses[[s[i]]] =  apply(x, 1, function(y) plotrix::std.error(y))
  }
  
  new_pdf_means = do.call(rbind, new_pdf_means)
  new_pdf_ses = do.call(rbind, new_pdf_ses)
  
  new_pdf_means = new_pdf_means[match(tree$tip.label, rownames(new_pdf_means)),]
  new_pdf_ses = new_pdf_ses[match(tree$tip.label, rownames(new_pdf_ses)),]
  
  #Calculate phylogenetic signal factoring in standard error
  out = list()
  print('calculating phylogenetic signal')
  for(i in 1:ncol(new_pdf_means)){
    #print(i)
    #out[[i]] = phylosig(tree, pdf_means[,i], se = pdf_ses[,i], nsim = 10000, test = TRUE)
    
    if(sum(new_pdf_means[,i]<quantile(unlist(new_pdf_means), prob = 0.3))>10){
      out[[i]] = list(K = NA)
    }else{
      out[[i]] = phylosig(tree, new_pdf_means[,i], se = new_pdf_ses[,i])
      #print(out[[i]]$K)
    }
    #print(out[[i]]$P)
  }
  
  l = list(new_pdf_means, as.numeric(na.omit(unlist(lapply(out, function(x) x$K)))))
  names(l) = c('pdf_means', 'ks')
  ks[[as.character(h)]] = l
  print(mean(ks[[as.character(h)]]$ks))
}

#Plot
plot(unlist(lapply(ks, function(x) mean(x$k))), 
     type = 'l', 
     ylim = c(0,1),
     ylab = 'Phylogenetic signal',
     xlab = 'PDF resolution (# of gridpoints)',
     bty = 'n',
     cex.axis = 1.5,
     cex.lab = 1.5,
     lwd = 1.5)

#Plot increasing resolution
par(mfrow = c(10,10), mar = c(0.25,0.25,0.25,0.25))
for(i in 1:length(ks)){
  image(matrix(ks[[i]]$pdf_means[2,], 
               ncol = as.numeric(names(ks)[i]), 
               as.numeric(names(ks)[i])),
        col = colorRampPalette(c('white', 'black'))(100),
        bty = 'n',
        xaxt = 'n',
        yaxt = 'n')
}

#######################################################
#####Compare distributions of phylogenetic signal######
#######################################################
#Combine into list
k = list(unlist(binary_k),
         as.numeric(na.omit(unlist(as.data.frame(occupancy_k)))),
         markov_k)
names(k) = c('binary', 'occupancy', 'markov')

x = c(jitter(rep(1, length(k[[1]]))),
      jitter(rep(2, length(k[[2]]))),
      jitter(rep(3, length(k[[3]]))))
y = unlist(k)

#Violin plot
par(bty = 'n')
vioplot::vioplot(k[[1]], 
                 k[[2]],
                 k[[3]],
                 cex.axis = 1.5,
                 cex.lab = 1.5,
                 cex.sub = 1.5,
                 ylab = 'Phylogenetic signal (K)',
                 las = 2,
                 names = c('Binary', 'Frequency', 'Transitions'),
                 col = c('cadetblue3', 'lightgoldenrod4', 'indianred3'),
                 border = c('cadetblue3', 'lightgoldenrod4', 'indianred3'))
abline(h = 1,
       lty = 'dashed',
       col = 'grey60')
title(main = paste('p =', signif(kruskal.test(k)$p.value, digits = 3)),
      cex.main = 1.5,
      font.main = 1.5)
cols = c('cadetblue3', 'lightgoldenrod4', 'indianred3')
for(i in c(1,2,3)){
  points(jitter(rep(i, length(k[[i]])), 1.5),
         k[[i]],
         pch = 20,
         col = darken_color(cols[i]))}

#Post hoc test
x = unlist(k)
g = c(rep('binary', length(k[[1]])),
      rep('freq', length(k[[2]])),
      rep('markov', length(k[[3]])))
d = dunn.test::dunn.test(x, g)

#Gap plot
par(bty = 'n')
gap.plot(x,
         y, 
         gap =c(1.5,6.5),
         ylim = c(0,7),
         ylab = 'Phylogenetic signal (K)',
         xlab = '',
         cex.axis = 1.5,
         cex.lab = 1.5,
         xaxt = 'n',
         pch = 20,
         cex = 1.5,
         col = alpha('grey70', 0.5),
         xlim = c(1, 3.5))
title(main = paste('p =', signif(kruskal.test(k)$p.value, digits = 3)),
      cex.main = 1.5,
      font.main = 1.5)
for(i in 1:length(k)){
  z = k[[i]]
  s = boxplot.stats(z)$stats
  points(i+0.5, s[3], cex = 3, pch = 20, col = 'grey50')
  segments(i+0.5, s[2], i+0.5, s[4], col = 'grey50', lwd = 5)
  segments(i+0.5, s[1], i+0.5, s[5], col = 'grey50', lwd = 2)
}

abline(h=seq(1.49,1.59,.001), col="white")
axis.break(2,1.5,style="slash", brw = 0.07)   
axis.break(4,1.5,style="slash", brw = 0.07) 

######################
#####Cophyloplots#####
######################
#Binary
d = dist(binary)
hcl = as.phylo(hclust(d))
obj = cophylo(tree, hcl)
plot(obj)

#Frequency
#pca = prcomp(pdf_means)
d = dist(pdf_means)
hcl = as.phylo(hclust(d))
obj = cophylo(tree, hcl)
plot(obj,
     link.type="curved",
     link.lwd=2,
     link.col=make.transparent("grey",0.7))

#Transitions
#pca = prcomp(trait)
d = dist(trait)
hcl = as.phylo(hclust(d))
obj = cophylo(tree, hcl)
plot(obj,
     link.type="curved",
     link.lwd=2,
     link.col=make.transparent("grey",0.7))

###############################
#####Bayestraits (via btw)#####
###############################
#Set working directory
setwd('bayes_traits/')

######Markov#####
#Create trait dataframe
pca = prcomp(trait)
dat = data.frame(species = rownames(trait),
                 trait = pca$x[,1:4])

#Run
setwd('bayes_traits/fossil_tree/markov_pc1/mod1/')
command_vec1 <- c("7", "2", 'VarRates', 'iterations 10010000')
mod1 = bayestraits(dat, 
                   tree, 
                   command_vec1,
                   remove_files = FALSE)

setwd('bayes_traits/fossil_tree/markov_pc1/mod2/')
command_vec2 <- c("7", "2", 'VarRates', 'iterations 10010000')
mod2 = bayestraits(dat, 
                   tree, 
                   command_vec2,
                   remove_files = FALSE)

setwd('bayes_traits/fossil_tree/markov_pc1/mod3/')
command_vec3 <- c("7", "2", 'VarRates', 'iterations 10010000')
mod3 = bayestraits(dat, 
                   tree, 
                   command_vec3,
                   remove_files = FALSE)

#Check for convergence
res1 <- mcmc(mod1$Log$results[,c(-1,-4)],
             start=min(mod1$Log$results$Iteration),
             end=max(mod1$Log$results$Iteration),thin=1000)

res2 <- mcmc(mod2$Log$results[,c(-1,-4)],
             start=min(mod2$Log$results$Iteration),
             end=max(mod2$Log$results$Iteration),thin=1000)

res3 <- mcmc(mod3$Log$results[,c(-1,-4)],
             start=min(mod3$Log$results$Iteration),
             end=max(mod3$Log$results$Iteration),thin=1000)

#Combine the three chains
res <- mcmc.list(res1,res2,res3)

#Trace and acf plots
traceplot(res[,c(1,3)])
acfplot(res[,c(1,3)])

#Get effective size (ideally > 200)
effectiveSize(res)

#Gelman-Rubin (should be ~1; if >1.05 run for longer)
gelman.diag(res,autoburnin=FALSE,multivariate=FALSE)

#Density plots of parameters across runs
densityplot(res[,-2])

#Choose best model using bayes factor
bftest(mod1, mod2)
bftest(mod1, mod3)
bftest(mod2, mod3)

#Process varrates output for best model
setwd('bayes_traits/fossil_tree/markov_pc1/mod1/')
markov_pc1_vrates_summary <- rjpp(logfile = 'data.txt.Log.txt',
                                  rjlog = 'data.txt.VarRates.txt', 
                                  rjtrees = 'data.txt.Output.trees',
                                  tree = tree)

#Plot
plot(markov_pc1_vrates_summary)

#####Occupancy by bin#####
#Create trait dataframe
pca = prcomp(pdf_means)
dat = data.frame(species = rownames(trait),
                 trait = pca$x[,1:13])

#Run
setwd('bayes_traits/fossil_tree/occupancy_pc1/mod1/')
command_vec1 <- c("7", "2", 'VarRates', 'iterations 50010000')
mod1 = bayestraits(dat, 
                   tree, 
                   command_vec1,
                   remove_files = FALSE)

setwd('bayes_traits/fossil_tree/occupancy_pc1/mod2/')
command_vec2 <- c("7", "2", 'VarRates', 'iterations 50010000')
mod2 = bayestraits(dat, 
                   tree, 
                   command_vec2,
                   remove_files = FALSE)

setwd('bayes_traits/fossil_tree/occupancy_pc1/mod3/')
command_vec3 <- c("7", "2", 'VarRates', 'iterations 50010000')
mod3 = bayestraits(dat, 
                   tree, 
                   command_vec3,
                   remove_files = FALSE)

#Check for convergence
res1 <- mcmc(mod1$Log$results[,c(-1,-4)],
             start=min(mod1$Log$results$Iteration),
             end=max(mod1$Log$results$Iteration),thin=1000)

res2 <- mcmc(mod2$Log$results[,c(-1,-4)],
             start=min(mod2$Log$results$Iteration),
             end=max(mod2$Log$results$Iteration),thin=1000)

res3 <- mcmc(mod3$Log$results[,c(-1,-4)],
             start=min(mod3$Log$results$Iteration),
             end=max(mod3$Log$results$Iteration),thin=1000)

#Combine the three chains
res <- mcmc.list(res1,res2,res3)

#Trace and acf plots
traceplot(res[,c(1,3)])
acfplot(res[,c(1,3)])

#Get effective size (ideally > 200)
effectiveSize(res)

#Gelman-Rubin (should be ~1; if >1.05 run for longer)
gelman.diag(res,autoburnin=FALSE,multivariate=FALSE)

#Density plots of parameters across runs
densityplot(res[,-2])

#Choose best model using bayes factor
bftest(mod1, mod2)
bftest(mod1, mod3)
bftest(mod2, mod3)

#Process varrates output for best model
setwd('bayes_traits/fossil_tree/occupancy_pc1/mod1/')
occupancy_pc1_vrates_summary <- rjpp(logfile = 'data.txt.Log.txt',
                                     rjlog = 'data.txt.VarRates.txt', 
                                     rjtrees = 'data.txt.Output.trees',
                                     tree = tree)

#Plot
plot(occupancy_pc1_vrates_summary)

#####Binary occupancy#####
#Create trait dataframe
pca = prcomp(binary)
dat = data.frame(species = rownames(trait),
                 trait = pca$x[,1:2])

#xy = unique(layout$xy_new)
#dat = lapply(strains, function(x) length(unique(x$xy_new))/length(xy))
#dat = dat[match(tree$tip.label, names(dat))]
#dat = data.frame(species = rownames(trait),
#trait = dat)

#Run
setwd('bayes_traits/fossil_tree/binary_pc1/mod1/')
command_vec1 <- c("7", "2", 'VarRates', 'iterations 10010000')
mod1 = bayestraits(dat, 
                   tree, 
                   command_vec1,
                   remove_files = FALSE)

setwd('bayes_traits/fossil_tree/binary_pc1/mod2/')
command_vec2 <- c("7", "2", 'VarRates', 'iterations 10010000')
mod2 = bayestraits(dat, 
                   tree, 
                   command_vec2,
                   remove_files = FALSE)

setwd('bayes_traits/fossil_tree/binary_pc1/mod3/')
command_vec3 <- c("7", "2", 'VarRates', 'iterations 10010000')
mod3 = bayestraits(dat, 
                   tree, 
                   command_vec3,
                   remove_files = FALSE)

#Check for convergence
res1 <- mcmc(mod1$Log$results[,c(-1,-4)],
             start=min(mod1$Log$results$Iteration),
             end=max(mod1$Log$results$Iteration),thin=1000)

res2 <- mcmc(mod2$Log$results[,c(-1,-4)],
             start=min(mod2$Log$results$Iteration),
             end=max(mod2$Log$results$Iteration),thin=1000)

res3 <- mcmc(mod3$Log$results[,c(-1,-4)],
             start=min(mod3$Log$results$Iteration),
             end=max(mod3$Log$results$Iteration),thin=1000)

#Combine the three chains
res <- mcmc.list(res1,res2,res3)

#Trace and acf plots
traceplot(res[,c(1,3)])
acfplot(res[,c(1,3)])

#Get effective size (ideally > 200)
effectiveSize(res)

#Gelman-Rubin (should be ~1; if >1.05 run for longer)
gelman.diag(res,autoburnin=FALSE,multivariate=FALSE)

#Density plots of parameters across runs
densityplot(res[,-2])

#Choose best model using bayes factor
bftest(mod1, mod2)
bftest(mod1, mod3)
bftest(mod2, mod3)

#Process varrates output for best model
setwd('bayes_traits/fossil_tree/binary_pc1/mod1/')
binary_vrates_summary <- rjpp(logfile = 'data.txt.Log.txt',
                              rjlog = 'data.txt.VarRates.txt', 
                              rjtrees = 'data.txt.Output.trees',
                              tree = tree)

#Plot
plotShifts(binary_vrates_summary)

#######################
#####Compare rates#####
#######################
#Calculate rates using Ronco et al. method
setwd('bayes_traits/fossil_tree/')
files = list.files()
files = files[grep('pc1', files)]

relRates = list()
for(i in 1:length(files)){
  print(i)
  setwd(paste(files[i], '/mod1', sep = ''))
  
  trees1 = read.nexus('data.txt.Output.trees')
  origtree = mytimetree
  
  out = matrix( NA, length(origtree$edge.length),length(trees1) )
  for (j in 1: dim(out)[2]) {
    out[,j] = trees1[[j]]$edge.length/ origtree$edge.length
  }
  
  r = rowMeans(out)
  
  relRates[[files[i]]] = r
  setwd('bayes_traits/fossil_tree/')
}

#Calculate mean trees 
#Mean trees
#setwd('bayes_traits/')
setwd('bayes_traits/fossil_tree/')
files = list.files()
files = files[grep('pc1', files)]
mean_trees = list()

for(i in 1:length(files)){
  setwd(paste(files[i], '/mod1', sep = ''))
  
  #####Calculate morphospace expansion/contraction#####
  #Get mean tree from Bayestraits (to run ancestral state reconstruction on)
  post_trees <- summariseTrees(reftree = tree,
                               trees = 'data.txt.Output.trees')
  mean_trees[[files[i]]] = post_trees$tree_summaries$mean_tree
  #setwd('bayes_traits/')
  setwd('bayes_traits/fossil_tree/')
}

par(mfrow = c(1,3))
lapply(mean_trees, function(x) plot(x))

#Phylogenies colored with rate
par(mfrow = c(1,3))
for(i in 1:length(relRates)){
  #if(i > 1){
  #rate = relRates[[i]]*100
  rate = relRates[[i]]
  cols = colorRampPalette(c("grey","gold","tomato","red", "darkred", "#3D0404"))(length(seq(0, round(max(rate), 3), 0.001)))
  names(cols) = seq(0, round(max(rate), 3), 0.001)
  cols = cols[match(round(rate, 3), names(cols))]
  #}else{
  #  rate = relRates[[i]]
  #  cols = colorRampPalette(c("grey","gold","tomato","red", "darkred", "#3D0404"))(length(seq(0, round(max(rate), 2), 0.01)))
  #  names(cols) = seq(0, round(max(rate), 2), 0.01)
  #  cols = cols[match(round(rate, 2), names(cols))]
  #}
  
  plot.phylo(mytimetree,
             edge.width = 1.5,
             edge.color = cols,
             show.tip.label = FALSE)
  #axisPhylo(cex.axis = 1.5, xlab = 'Time (mya)')
  title(main = unlist(lapply(strsplit(names(relRates)[i], '_'), function(y) y[1])),
        cex.main = 1.5,
        font.main = 1)
  
  #tip.color = dros_cols[match(mytimetree$tip.label, names(dros_cols))])
}

names(relRates) = c('binary', 'markov', 'occupancy')

#Compare rates across traits
markov_rates = markov_pc1_vrates_summary$scalars$meanRate[-1]/max(markov_pc1_vrates_summary$scalars$meanRate[-1])
occupancy_rates = occupancy_pc1_vrates_summary$scalars$meanRate[-1]/max(occupancy_pc1_vrates_summary$scalars$meanRate[-1])
binary_rates = binary_vrates_summary$scalars$meanRate[-1]/max(binary_vrates_summary$scalars$meanRate[-1])

#Plot distributions
#Violin plot
par(bty = 'n')
vioplot::vioplot(relRates[[1]]/max(relRates[[1]]), 
                 relRates[[3]]/max(relRates[[3]]),
                 relRates[[2]]/max(relRates[[2]]),
                 cex.axis = 1.5,
                 cex.lab = 1.5,
                 cex.sub = 1.5,
                 ylab = 'Relative rate of evolution',
                 las = 2,
                 names = c('Binary', 'Frequency', 'Transitions'),
                 col = c('cadetblue3', 'lightgoldenrod4', 'indianred3'),
                 border = c('cadetblue3', 'lightgoldenrod4', 'indianred3'))
title(main = paste('p =', signif(kruskal.test(relRates)$p.value, digits = 3)),
      cex.main = 1.5,
      font.main = 1.5)

for(i in c(1,3,2)){
  points(jitter(rep(i, length(relRates[[i]])), 1.5),
         relRates[[i]]/max(relRates[[i]]),
         pch = 20,
         col = 'grey80')}

#Dunn's test
dunn.test::dunn.test(relRates)

#Beeswarm
par(mfrow = c(1,3))
beeswarm(unlist(lapply(relRates, function(x) x/max(x)))~c(rep(names(relRates[1]), length(relRates[[1]])),
                                                          rep(names(relRates[2]), length(relRates[[2]])),
                                                          rep(names(relRates[3]), length(relRates[[3]]))),
         pch = 20,
         col = 'grey80',
         ylab = 'Relative rate of evolution',
         xlab = '',
         cex.lab = 1.5,
         cex.axis = 1.5,
         spacing = 0.75,
         las = 2,
         bty = 'n',
         method = 'center')
title(main = paste('p =', signif(kruskal.test(relRates)$p.value, digits = 3)),
      cex.main = 1.5,
      font.main = 1.5)
for(i in 1:length(relRates)){
  x = relRates[[i]]/max(relRates[[i]])
  s = boxplot.stats(x)$stats
  points(i+0.5, s[3], cex = 3, pch = 20, col = 'grey50')
  segments(i+0.5, s[2], i+0.5, s[4], col = 'grey50', lwd = 5)
  segments(i+0.5, s[1], i+0.5, s[5], col = 'grey50', lwd = 2)
}

###############################################################################################################################
#####Morphospace expansion/packing analyses
#####(following https://github.com/cichlidx/ronco_et_al/blob/master/trait_evolution/03_TraitEvolution/scripts/05_AncRec.R)#####
###############################################################################################################################
setwd('bayes_traits/fossil_tree/')
files = list.files()
morpho_results = list()

data = list()
data$binary = data.frame(species = rownames(binary),
                         trait = prcomp(binary)$x[,1:2])
data$markov = data.frame(species = rownames(trait),
                         trait = prcomp(trait)$x[,1:4])
data$occupancy = data.frame(species = rownames(pdf_means),
                            trait = prcomp(pdf_means)$x[,1:13])

for(h in 1:length(files)){
  
  #Set to directory with results
  setwd(paste(files[h], '/mod1', sep = ''))
  
  #Create empty list to save results
  l = list()
  
  #####Calculate morphospace expansion/contraction#####
  #Get mean tree from Bayestraits (to run ancestral state reconstruction on)
  post_trees <- summariseTrees(reftree = tree,
                               trees = 'data.txt.Output.trees')
  mean_tree = post_trees$tree_summaries$mean_tree
  
  #Run ml ancestral state reconstruction on PC1 using phytools
  x = data[[h]][,2]
  names(x) = data[[h]][,1]
  anc_tree = fastAnc(mean_tree, x)
  
  #Add ancestral states and tips to original tree
  new = c(data[[h]][,2], anc_tree)
  names(new) = c(data[[h]][,1], names(anc_tree))
  
  #Get times per nodes
  H = nodeHeights(mytimetree)
  root = max(H[,2])
  
  #Round root
  root2 = floor(root*1000)/1000
  
  #Define time slices
  STEP = 0.1
  n = STEP
  v = c(seq(from=0, to=root2, by= n) , root2)
  
  #Combine states of the nodes with nodeages
  xx = node.depth.edgelength(mytimetree)  
  yy.bayPC = cbind(new, xx)          
  
  #For each timepoint reconstruct the state linearly to the distance between the two node stages
  Yt = list()
  for ( i in c(1: (length(v)-1)) ) {      
    tmp.bayPC3      = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
    tmp.bayPC.node3 = mytimetree$edge[tmp.bayPC3,]                   ### get nodes of all branches between root and timepoint 
    start           = yy.bayPC[tmp.bayPC.node3[,1],]
    end             = yy.bayPC[tmp.bayPC.node3[,2],]
    t               = rep(v[i],times= length(end[,1]))
    y3              =(((end[,1]-start[,1])/(end[,2]-start[,2])) * (v[i] - start[,2])) + start[,1]
    Yt[[i]]         =cbind(y3, t)
  }
  
  #Clean
  v2 = v[-(length(v)-1)]
  
  out3=NULL
  for ( k in 1: length(v2)){
    tmp3= unlist(Yt[[k]])
    tmp3.2= data.frame(tmp3, "sample" = v2[k])
    out3= rbind(out3, tmp3.2)
  }
  
  #Simplify
  d=out3[,c(1,3)]
  
  #Convert sample to factor
  d$sample2 = as.factor(d$sample)
  
  #Make function for normalization
  normalize = function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
  
  #Create empty dataeframe to collect stats
  statsOtime = as.data.frame(matrix(NA, length(levels(d$sample2)) , 1 ))
  statsOtime[,1] = as.numeric(as.character(levels(d$sample2)))
  names(statsOtime) = "time"
  
  #Get range from each timeslice:
  minv = aggregate(d[,1], list(d$sample2), min)
  maxv = aggregate(d[,1], list(d$sample2), max)
  
  #Range of the axes
  statsOtime$range = abs(minv[,-1]-maxv[,-1])
  
  #Get number of linegaes per timeslice
  count_per_bin = aggregate(d[,1], list(d$sample2), length)
  
  #Combine and remove the first row = root value as diversity= 0
  statsOtime$numLineages = count_per_bin[,2]
  
  #Function to plot
  plot_morphospace_expansions = function(slices, root, step = 0.15){
    
    nbins = trunc(root/step)
    
    a = hist2d(slices$sample, slices$y3, nbins= nbins, axes=F, show=F)
    maxa = max(a$counts)
    tck = -0.02
    xlim=c(0-(root2-root),root2)
    
    hm_col_scale = colorRampPalette(c( "grey50", "black"))(maxa)
    hm_col_scale2= c("white", hm_col_scale)
    hist2d(slices$sample, slices$y3, nbins= nbins, col= hm_col_scale2, axes=F, xlim=xlim)
    axis(1, at = rev(root2- c(seq(0, root2, 1)))-(root2-root), rev(c(seq(0, root2, 1))), 
         tck = tck, 
         las = 1, 
         cex.axis = 1.5)
    axis(2, 
         tck = tck, 
         las = 1, 
         cex.axis = 1.5)
    mtext("Time (Ma)", 1, line=1.4, cex=1.5, padj = 1.5)
  }
  
  #Plot
  plot_morphospace_expansions(d, root)
  
  #Add marker of where 50% of species have evolved
  abline(v = 37, lwd = 5, col = alpha('cyan4', 0.25))
  
  l$morphospace_slices = d
  l$statsoftime = statsOtime
  
  #####Calculate mean rate per slice#####
  #get trees
  trees1 = read.nexus('data.txt.Output.trees')
  
  #Ladderize the trees
  origtree = ladderize(mytimetree)
  trees1 = lapply(trees1, ladderize)
  class(trees1) = "multiPhylo"
  
  #Get relative rates
  out= matrix( NA, length(origtree$edge.length),length(trees1) )
  for ( i in 1: dim(out)[2]) {
    out[,i] = trees1[[i]]$edge.length/origtree$edge.length
  }
  
  Yt = matrix( NA,length(v) ,dim(out)[2])
  for ( i in c(1: length(v)) ) {
    tmp.bayPC3 = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
    tmp_rates= out[tmp.bayPC3,]     ###  find those edges in the tree data:
    if (length(tmp.bayPC3)>=2) {
      Yt[i,] = apply(tmp_rates, 2, mean)
    } else {
      Yt[i,] = mean(tmp_rates)        
    }
  }
  
  #Extract mean and standard deviation
  mean_rate= rowMeans(Yt, na.rm=T)
  tmp = apply(Yt ,1, sd)
  yy = c(mean_rate+tmp, rev(mean_rate-tmp))
  xx = c(v , rev(v))
  polySD = cbind(x=xx, y=yy)
  
  #Clean
  mean_rates= cbind(v,mean_rate)
  
  #Convert into objects for plot
  mean_rateTT = data.frame(v = mean_rates[,1],
                           mean_rate = mean_rates[,2])
  CI_rateTT = polySD
  
  #Plot
  xlim=c(0-(root2-root),root2)
  ylim= c(0, max(CI_rateTT[, "y"]) *1.25)
  tck = -0.02
  
  plot(mean_rateTT$v, mean_rateTT$mean_rate, ylim=ylim,xlim=xlim ,t="n", axes=F, xlab="", ylab="")
  polygon(CI_rateTT[, "x"], CI_rateTT[, "y"], col = alpha("deepskyblue4",0.2), border = NA)
  lines(mean_rateTT$v, mean_rateTT$mean_rate, type="l",col="deepskyblue4", lwd=3)
  axis(1,at= rev(root2- c(seq(0, root2, 1)))-(root2-root), rev(c(seq(0, root2, 1))), tck= tck, cex.axis= 1.5)
  axis(2, tck= tck, las=1, cex.axis= 1.5)
  mtext("Time (Ma)", 1, line=1.4, cex=1.5, padj = 1.5)
  mtext("Relative rate of evolution", 2, line=1, cex=1.5)
  abline(h=1, col="grey", lty="dashed")
  
  l$rate_mean = mean_rateTT
  l$rate_ci = CI_rateTT
  
  morpho_results[[files[h]]] = l
  setwd('bayes_traits/fossil_tree/')
}

#Save
saveRDS(morpho_results, 'binary_markov_occup_morphospace_expansion_results_120720.RDS')
saveRDS(morpho_results, 'binary_markov_occup_morphospace_expansion_results_fossiltree_050321.RDS')

#Load
morpho_results = readRDS('binary_markov_occup_morphospace_expansion_results_120720.RDS')
morpho_results = readRDS('binary_markov_occup_morphospace_expansion_results_fossiltree_050321.RDS')

#Plot relative rate over time (with sd as in Ronco)
par(mfrow = c(1,3))
for(i in c(1,3,2)){
  CI_rateTT = morpho_results[[i]]$rate_ci
  mean_rateTT = morpho_results[[i]]$rate_mean
  
  #CI_rateTT[,2] = CI_rateTT[,2]/max(mean_rateTT[,2])
  #mean_rateTT[,2] = mean_rateTT[,2]/max(mean_rateTT[,2])
  
  ylim= c(0, max(CI_rateTT[, "y"]) *1.25)
  
  #Get times per nodes
  H = nodeHeights(mytimetree)
  root = max(H[,2])
  
  #Round root
  root2 = floor(root*1000)/1000
  
  xlim=c(0-(root2-root),root2)
  #ylim = c(min(CI_rateTT[,2]), max(CI_rateTT[,2]))
  #ylim = c(0,1)
  tck = -0.02
  
  plot(mean_rateTT$v, mean_rateTT$mean_rate, ylim=ylim,xlim=xlim ,t="n", axes=F, xlab="", ylab="")
  polygon(CI_rateTT[, "x"], CI_rateTT[, "y"], col = alpha("deepskyblue4",0.2), border = NA)
  lines(mean_rateTT$v, mean_rateTT$mean_rate, type="l",col="deepskyblue4", lwd=3)
  axis(1,at= rev(root2- c(seq(0, root2, 1)))-(root2-root), rev(c(seq(0, root2, 1))), tck= tck, cex.axis= 1.5)
  axis(2, tck= tck, las=1, cex.axis= 1.5)
  mtext("Time (Ma)", 1, line=1.4, cex=1.5, padj = 1.5)
  mtext("Relative rate of evolution", 2, line=1, cex=1.5)
  abline(h=1, col="grey", lty="dashed")
}

#Plot just normalized mean
plot(morpho_results$binary_pc1$rate_mean$mean_rate/max(morpho_results$binary_pc1$rate_mean$mean_rate), 
     type = 'l',
     ylim = c(0,1),
     ylab = 'Relative rate of evolution',
     cex.axis = 1.5,
     cex.lab = 1.5)
lines(morpho_results$occupancy_pc1$rate_mean$mean_rate/max(morpho_results$occupancy_pc1$rate_mean$mean_rate), 
      col = 'orange')
lines(morpho_results$markov_pc1$rate_mean$mean_rate/max(morpho_results$markov_pc1$rate_mean$mean_rate),
      col = 'green')

#Plot species accumulation
plot(morpho_results$binary_pc1$statsoftime$time-40, 
     morpho_results$binary_pc1$statsoftime$numLineages,
     type = 'l',
     lwd = 1.5,
     cex.lab = 1.5,
     cex.axis = 1.5,
     xlab = 'Time (mya)',
     ylab = '# species/strains')

#Morphospace filling in 2d
setwd('bayes_traits/fossil_tree/')
files = list.files()
morpho_2d_filling = list()
for(h in 1:length(files)){
  
  print(files[h])
  #Set to directory with results
  setwd(paste(files[h], '/mod1', sep = ''))
  
  #Create empty list to save results
  l = list()
  
  #####Calculate morphospace expansion/contraction#####
  #Get mean tree from Bayestraits (to run ancestral state reconstruction on)
  post_trees <- summariseTrees(reftree = tree,
                               trees = 'data.txt.Output.trees')
  mean_tree = post_trees$tree_summaries$mean_tree
  
  pca_ancs = list()
  for(b in 1:2){
    
    #Run ml ancestral state reconstruction on PC1 using phytools
    x = data[[h]][,(b+1)]
    names(x) = data[[h]][,1]
    anc_tree = fastAnc(mean_tree, x)
    
    #Add ancestral states and tips to original tree
    new = c(data[[h]][,(b+1)], anc_tree)
    names(new) = c(data[[h]][,1], names(anc_tree))
    
    #Get times per nodes
    H = nodeHeights(mytimetree)
    root = max(H[,2])
    
    #Round root
    root2= floor(root*1000)/1000
    
    #Define time slices
    STEP = 0.1
    n = STEP
    v = c(seq(from=0, to=root2, by= n) , root2)
    
    #Combine states of the nodes with nodeages
    xx = node.depth.edgelength(mytimetree)  
    yy.bayPC = cbind(new, xx)          
    
    #For each timepoint reconstruct the state linearly to the distance between the two node stages
    Yt = list()
    for ( i in c(1: (length(v)-1)) ) {      
      tmp.bayPC3      = intersect(which(H[,1]<= v[i] ), which( H[,2] >= v[i] ))  ###  get row number for the nodesHeights within root and timepoint  ---   get all edges crossing the timeslcie v:tips, but not ending before ( not adding nodes which are not present anymore)
      tmp.bayPC.node3 = mytimetree$edge[tmp.bayPC3,]                   ### get nodes of all branches between root and timepoint 
      start           = yy.bayPC[tmp.bayPC.node3[,1],]
      end             = yy.bayPC[tmp.bayPC.node3[,2],]
      t               = rep(v[i],times= length(end[,1]))
      y3              =(((end[,1]-start[,1])/(end[,2]-start[,2])) * (v[i] - start[,2])) + start[,1]
      Yt[[i]]         =cbind(y3, t)
    }
    
    #Clean
    v2 = v[-(length(v)-1)]
    
    out3=NULL
    for ( k in 1: length(v2)){
      tmp3= unlist(Yt[[k]])
      tmp3.2= data.frame(tmp3, "sample" = v2[k])
      out3= rbind(out3, tmp3.2)
    }
    
    #Simplify
    d=out3[,c(1,3)]
    
    #Convert sample to factor
    d$sample2 = as.factor(d$sample)
    
    #Make function for normalization
    normalize = function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T))}
    
    #Create empty dataeframe to collect stats
    statsOtime = as.data.frame(matrix(NA, length(levels(d$sample2)) , 1 ))
    statsOtime[,1] = as.numeric(as.character(levels(d$sample2)))
    names(statsOtime) = "time"
    
    #Get range from each timeslice:
    minv = aggregate(d[,1], list(d$sample2), min)
    maxv = aggregate(d[,1], list(d$sample2), max)
    
    #Range of the axes
    statsOtime$range = abs(minv[,-1]-maxv[,-1])
    
    #Get number of linegaes per timeslice
    count_per_bin = aggregate(d[,1], list(d$sample2), length)
    
    #Combine and remove the first row = root value as diversity= 0
    statsOtime$numLineages = count_per_bin[,2]
    
    #Plot
    plot_morphospace_expansions(d, root)
    
    #Add marker of where 50% of species have evolved
    abline(v = 37, lwd = 5, col = alpha('cyan4', 0.25))
    
    l$morphospace_slices = d
    l$statsoftime = statsOtime
    
    pca_ancs[[b]] = d
  }
  morpho_2d_filling[[files[h]]] = pca_ancs
  setwd('bayes_traits/fossil_tree/')
}

#Plot
#Get times per nodes
H = nodeHeights(mytimetree)
root = max(H[,2])

#Round root
root2 = floor(root*1000)/100

par(mfrow = c(4,3))
#Phylogenies colored with rate
for(i in 1:length(relRates)){
  #if(i > 1){
  #rate = relRates[[i]]*100
  rate = relRates[[i]]
  cols = colorRampPalette(c("grey","gold","tomato","red", "darkred", "#3D0404"))(length(seq(0, round(max(rate), 3), 0.001)))
  names(cols) = seq(0, round(max(rate), 3), 0.001)
  cols = cols[match(round(rate, 3), names(cols))]
  
  plot.phylo(mytimetree,
             edge.width = 1.5,
             edge.color = cols,
             show.tip.label = FALSE)
  title(main = unlist(lapply(strsplit(names(relRates)[i], '_'), function(y) y[1])),
        cex.main = 1.5,
        font.main = 1)
  
}
for(i in 1:length(morpho_2d_filling)){
  if(i==2){
    pcs = plot_morpho_filling_2d(pc1 = round(morpho_2d_filling[[i]][[1]]$y3, 3),
                                 pc2 = round(morpho_2d_filling[[i]][[2]]$y3, 3),
                                 dates = round(morpho_2d_filling[[i]][[1]]$sample, 2),
                                 #contour = TRUE,
                                 normalize_pcs = TRUE)
  }else{
    pcs = plot_morpho_filling_2d(pc1 = round(morpho_2d_filling[[i]][[1]]$y3, 1),
                                 pc2 = round(morpho_2d_filling[[i]][[2]]$y3, 1),
                                 dates = round(morpho_2d_filling[[i]][[1]]$sample, 2),
                                 #contour = TRUE,
                                 normalize_pcs = TRUE)
  }
  
  #cols = dros_cols[match(rownames(data$binary), names(dros_cols))]
  #points(data[[i]][,2:3], pch = 21, bg = cols, col = 'black', cex = 1.5)
  
  title(main = unlist(lapply(strsplit(names(morpho_2d_filling)[i], '_'), function(y) y[1])),
        cex.main = 1.5,
        font.main = 1)
}
fill = list()
for(i in 1:length(morpho_2d_filling)){
  pcs = plot_morpho_filling_2d(pc1 = round(morpho_2d_filling[[i]][[1]]$y3, 2),
                               pc2 = round(morpho_2d_filling[[i]][[2]]$y3, 2),
                               dates = round(morpho_2d_filling[[i]][[1]]$sample, 2),
                               plot = FALSE)
  pcs$pc1 = pcs$pc1+(abs(min(pcs$pc1)))
  pcs$pc1 = pcs$pc1/max(pcs$pc1)
  pcs$pc2 = pcs$pc2+(abs(min(pcs$pc2)))
  pcs$pc2 = pcs$pc2/max(pcs$pc2)
  
  #date = round(pcs$date, 1)
  #time = seq(0, 41, 0.1)
  time = seq(0, 47.75, 0.01)
  #d = cumsum(date)/max(cumsum(date))
  #d = d[match(time, date)]
  #names(d) = time
  #d = zoo::na.locf(d)
  
  cov = list()
  for(j in 1:nrow(pcs)){
    tmp = pcs[1:j,]
    cov[[j]] = compute.coverage(min(pcs$pc1), max(pcs$pc1), min(pcs$pc2), max(pcs$pc2),
                                boxsize = 0.1, x = tmp$pc1, y = tmp$pc2)
  }
  cov = unlist(lapply(cov, function(x) x[3]))
  names(cov) = pcs$date
  cov = cov[!duplicated(names(cov))]
  d = cov[match(time, names(cov))]
  names(d) = time
  d = zoo::na.locf(d)
  
  fill[[i]] = d
  
  cols = colorRampPalette(c(colorRampPalette(c('black', 'midnightblue'))(50), pals::parula(25)))(length(seq(0, 47.75, 0.01)))
  names(cols) = seq(0, 47.75, 0.01)
  cols = cols[match(time, names(cols))]
  
  plot(d,
       #type = 'n',
       cex.axis = 1.5,
       type = 'h',
       col = cols,
       cex.lab = 1.5,
       xlab = 'Time (mya)',
       ylab = '% of morphospace filled',
       xaxt = 'n',
       ylim = c(0, 0.4))
  for(j in 1:length(d)){
    segments(j, d[j], j+1, d[j+1],
             col = cols[j],
             lwd = 3)}
  axis(1, 
       at = seq(0, 4776, 477.6), 
       labels = rev(seq(0, 41, 4)),
       las = 1, 
       cex.axis = 1.5)
}

#Plot PCs as morphospace density over time
#Get amount of space filled
#par(mfrow = c(1,3))
for(h in 1:3){
  #x = morpho_2d_filling[[h]][[1]]
  #x = split(x, x$sample)
  #d = c()
  #for(i in 1:length(x)){
  #  d = c(d, abs(max(x[[i]]$y3)-min(x[[i]]$y3)))
  #}
  d = fill[[h]]
  d = d/max(d)
  d1 = as.numeric(names(d)[which.min(abs(d - 0.5))])
  d2 = as.numeric(names(d)[which.min(abs(d - 0.25))])
  
  cols = colorRampPalette(c(colorRampPalette(c('black', 'midnightblue'))(50), pals::parula(25)))(length(seq(0, max(morpho_2d_filling[[h]][[1]]$sample), 0.1)))
  names(cols) = seq(0, max(morpho_2d_filling[[h]][[1]]$sample), 0.1)
  cols = cols[match(time, names(cols))]
  cols = cols[match(morpho_2d_filling[[h]][[1]]$sample, names(cols))]
  
  plot(morpho_2d_filling[[h]][[1]]$sample,
       morpho_2d_filling[[h]][[1]]$y3, 
       pch = 20, 
       #col = alpha('grey40', 0.5),
       col = cols,
       ylab = 'PC1',
       xlab = 'Time (Ma)',
       cex.axis = 1.5,
       cex.lab = 1.5,
       xaxt = 'n',
       bty = 'n')
  abline(v = d1, col = 'red')
  abline(v = d2, col = 'red')
  
  axis(1, 
       at = seq(0, 47.8, 4.78), 
       labels = rev(seq(0, 41, 4)),
       las = 1, 
       cex.axis = 1.5)
  
  title(main = paste(names(morpho_2d_filling)[h], 
                     round(max(morpho_2d_filling[[h]][[1]]$sample), 2)-d1), 
        cex.main = 1.5, 
        font.main = 1)
}

#Legend
cols = colorRampPalette(col)(length(seq(0, max(pcs$date), 0.01)))
names(cols) = seq(0, max(pcs$date), 0.01)
plot(seq(1,length(cols),1),
     rep(1, length(cols)),
     col = NULL,
     ylim = c(0,1),
     xaxt = 'n',
     yaxt = 'n',
     ylab = '',
     xlab = '',
     cex.lab = 1.5,
     bty = 'n')
for(i in 1:(length(cols)-1)){
  rect(i, 0, i+1, 1, col = cols[i], border = cols[i])
}
axis(1, 
     at = seq(0, 60, 10), 
     labels = seq(0, 6, 1),
     cex.axis = 1.5)
