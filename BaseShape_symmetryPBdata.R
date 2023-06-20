source("E:/ProBound/PlotLogo_function.R")
if (!require('RJSONIO')){
  install.packages('RJSONIO')
  require('RJSONIO')
}
if (!require('plotly')){
  install.packages('plotly')
  require('plotly')
}
if (!require("processx")){
  install.packages("processx")
  require("processx")
}
if (!require("umap")){
  install.packages("umap")
  require("umap")
}
source('E:/Desktop/HB_Rotation_1/Family_code/FCfunctions.R')
library("rjson")
library(proteinStructureBoost)

####load mono-di matrix for probound_run results ####
proBound_runDIR <- 'E:\\Desktop\\HB_Rotation_1\\ProBound_run\\'
proBound_files <- list.files(proBound_runDIR)
modelFile_Template <- 'E:/Desktop/HB_Rotation_1/ProBound_run/$modelFile$/result/fit.models.consensus.json'
rows <- c("A" ,   "C"  ,  "G"   , "T"  ,  "AA:1", "AC:1" ,"AG:1", "AT:1", "CA:1", "CC:1",
          "CG:1", "CT:1" ,"GA:1", "GC:1", "GG:1", "GT:1", "TA:1", "TC:1", "TG:1", "TT:1")
rec_seq <- 'CANNTG'
blank <- data.frame(P00 = c(1,1,1,1), P0 = c(1,1,1,1))
di.blank <- rep(0,16)
monodi_motifs <- list()
####get quality####
bHLH_info <- read.csv(file = "E:/Desktop/HB_Rotation_1/Base_Shape_rec/bHLH_info.csv", row.names = 1)
bHLH_index <- bHLH_info[,1:2]
bHLH_index$quality <- 0
for(i in 1:nrow(bHLH_index)){
  tryCatch({
    dir <- paste0(proBound_runDIR,bHLH_index$gene_symbol[i],'_',bHLH_index$study[i])
    quality <- readChar(paste0(dir,'//fitEval.txt'),nchars = 100)
    quality <- strsplit(quality, ', ')[[1]][1]
    quality <- strsplit(quality, ' = ')[[1]][2]
    bHLH_index$quality[i] <- quality
  }, error = function(e){})
}
ProBound_quality <- bHLH_index
ProBound_quality <- ProBound_quality[ProBound_quality$study != 'Isakova2017',]
n <- 0
bHLH_index <- bHLH_index[bHLH_index$study != 'Isakova2017',]
mono_motifs <- list()
motif_model <- c()
symmetry <- c()
model_score <- c()
for(i in 1:nrow(bHLH_index)){
  tryCatch({
    proBound_file <- paste0(bHLH_index$gene_symbol[i], '_',bHLH_index$study[i] )
    modelFile <- gsub("\\$modelFile\\$", proBound_file, modelFile_Template)
    name <- proBound_file
    result <- parse_json(file = modelFile)
    mono <- result[[2]]$mononucleotide
    mono.df <- cbind.data.frame(name = rows[1:4], mono)
    
    JSON_matrix <- exp(mono)
    #JSON_matrix <- cbind.data.frame(blank, JSON_matrix, blank)
    bind_pos <- find_binding_site(JSON_matrix, rec_seq)$max_pos[1]
    score1 <- max(find_binding_site(JSON_matrix, rec_seq)$scores)
    JSON_matrix_motif1 <- JSON_matrix[,c(bind_pos:(bind_pos+nchar(rec_seq)-1))]
    
    mono <- result[[3]]$mononucleotide
    if(length(mono) != 0){
      mono.df <- cbind.data.frame(name = rows[1:4], mono)
      
      JSON_matrix <- exp(mono)
      #JSON_matrix <- cbind.data.frame(blank, JSON_matrix, blank)
      bind_pos <- find_binding_site(JSON_matrix, rec_seq)$max_pos[1]
      score2 <- max(find_binding_site(JSON_matrix, rec_seq)$scores)
      JSON_matrix_motif2 <- JSON_matrix[,c(bind_pos:(bind_pos+nchar(rec_seq)-1))]
    }else{
      score2 <- -Inf
    }
    
    
    if(score2 <= score1){
      JSON_matrix_motif <- JSON_matrix_motif1
      motif_model <- c(motif_model, 1)
      model_score <- c(model_score, score1)
    }else{
      JSON_matrix_motif <- JSON_matrix_motif2
      motif_model <- c(motif_model, 2)
      model_score <- c(model_score, score2)
    }
    pos_index <- c('P-3','P-2','P-1','P1','P2','P3')
    colnames(JSON_matrix_motif) <- pos_index
    
    if(sum(JSON_matrix_motif[,3]) == sum(JSON_matrix_motif[,4])){
      symmetry <- c(symmetry, 1)
    }else{
      symmetry <- c(symmetry, 0)
    }
    
    n <- n + 1
    mono_motif <- list(name = name, monoMatrix = JSON_matrix_motif)
    mono_motifs[[n]] <- mono_motif
    print(n)
  }, error = function(e){})
}
bHLH_model_info <- data.frame(bHLH_index, motif_model, model_score, symmetry, ProBound_quality$quality)
bHLH_sym <- bHLH_model_info[bHLH_model_info$symmetry == 1,]
bHLH_names <- unique(bHLH_sym$gene_symbol)

bHLH_index <- data.frame(NULL)
for(i in 1:length(bHLH_names)){
  qualities <- bHLH_sym[bHLH_sym$gene_symbol == bHLH_names[i],]
  add <- qualities[qualities$quality == max(qualities$quality),]
  bHLH_index <- rbind.data.frame(bHLH_index, add[1,])
}
nrow(bHLH_index)

bHLH_index <- bHLH_index[bHLH_index$quality > 0.15,]
#reenter with refined standards
n <- 0
mono_motifs <- list()
for(i in 1:nrow(bHLH_index)){
  tryCatch({
    proBound_file <- paste0(bHLH_index$gene_symbol[i], '_',bHLH_index$study[i] )
    modelFile <- gsub("\\$modelFile\\$", proBound_file, modelFile_Template)
    name <- proBound_file
    result <- parse_json(file = modelFile)
    mono <- result[[2]]$mononucleotide
    mono.df <- cbind.data.frame(name = rows[1:4], mono)
    
    JSON_matrix <- exp(mono)
    #JSON_matrix <- cbind.data.frame(blank, JSON_matrix, blank)
    bind_pos <- find_binding_site(JSON_matrix, rec_seq)$max_pos[1]
    score1 <- max(find_binding_site(JSON_matrix, rec_seq)$scores)
    JSON_matrix_motif1 <- JSON_matrix[,c(bind_pos:(bind_pos+nchar(rec_seq)-1))]
    
    mono <- result[[3]]$mononucleotide
    if(length(mono) != 0){
      mono.df <- cbind.data.frame(name = rows[1:4], mono)
      
      JSON_matrix <- exp(mono)
      #JSON_matrix <- cbind.data.frame(blank, JSON_matrix, blank)
      bind_pos <- find_binding_site(JSON_matrix, rec_seq)$max_pos[1]
      score2 <- max(find_binding_site(JSON_matrix, rec_seq)$scores)
      JSON_matrix_motif2 <- JSON_matrix[,c(bind_pos:(bind_pos+nchar(rec_seq)-1))]
    }else{
      score2 <- -Inf
    }
    
    
    if(score2 <= score1){
      JSON_matrix_motif <- JSON_matrix_motif1
      motif_model <- c(motif_model, 1)
      model_score <- c(model_score, score1)
    }else{
      JSON_matrix_motif <- JSON_matrix_motif2
      motif_model <- c(motif_model, 2)
      model_score <- c(model_score, score2)
    }
    pos_index <- c('P-3','P-2','P-1','P1','P2','P3')
    colnames(JSON_matrix_motif) <- pos_index
    
    if(sum(JSON_matrix_motif[,3]) == sum(JSON_matrix_motif[,4])){
      symmetry <- c(symmetry, 1)
    }else{
      symmetry <- c(symmetry, 0)
    }
    
    n <- n + 1
    mono_motif <- list(name = name, monoMatrix = JSON_matrix_motif)
    mono_motifs[[n]] <- mono_motif
    print(n)
  }, error = function(e){})
}



bHLH_pb <- read.table("E:/Desktop/HB_Rotation_1/Base_Shape_rec/bHLH_pb.sto.stockholm", quote="\"", fill = T)
bHLH_pb <- bHLH_pb[-nrow(bHLH_pb),]
names <- unique(bHLH_pb$V1)
bHLH_pbAlignment <- data.frame(NULL)
for(i in 1:length(names)){
  ali <- bHLH_pb[bHLH_pb$V1 == names[i],2]
  alignment <- paste(ali,collapse = '')
  add <- data.frame(name = names[i], alignment = alignment)
  bHLH_pbAlignment <- rbind.data.frame(bHLH_pbAlignment, add)
}
bHLH_pbAlignment$alignment <- substr(bHLH_pbAlignment$alignment, 694, (694+55))



bHLH_motifs <- list()

for(i in 1:length(mono_motifs)){
  add <- list(name = mono_motifs[[i]]$name, matrix = mono_motifs[[i]]$monoMatrix)
  bHLH_motifs[[length(bHLH_motifs)+1]] <- add

}


pos_index <- c('P-3','P-2','P-1','P1','P2','P3')
pVal.Pos <- data.frame(NULL)
for(AApos in 1:nchar(bHLH_pbAlignment$alignment[1])){
  ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
  PosAA <- c()
  for(i in 1:length(bHLH_motifs)){
    PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
  }
  if(length(unique(PosAA)) == 1){
    addLine <- NA
    pVal.Pos <- rbind.data.frame(pVal.Pos, addLine)
    next()
  }
  addLine <- c()
  for(n.pos in 1:length(pos_index)){
    pos_matrix <- gene2pos(bHLH_motifs, pos = pos_index[n.pos])
    df <- matrix2tetrahedron(pos_matrix)
    
    df <- data.frame(AA = PosAA, A = df[,1], C = df[,2], G = df[,3])
    man.test <- manova(cbind(A,C,G) ~ AA, data = df)
    
    #test <- lm(A+G+C ~ AA, data = df)
    #man.test <- car::Anova(test, type = 'II')
    
    summary <- summary(man.test, tol=0)
    p.val <- summary$stats[1,6]
    #p.val <- man.test$`Pr(>F)`[2]
    addLine <- c(addLine, p.val)
  }
  pVal.Pos <- rbind.data.frame(pVal.Pos, addLine)
}
rownames(pVal.Pos) <- paste('AA', c(1:nchar(bHLH_pbAlignment$alignment[1])), sep = '')
colnames(pVal.Pos) <- pos_index
pVal.Pos <- -log(pVal.Pos,10)
#write.csv(pVal.Pos, file = 'E:/Desktop/HB_Rotation_1/Base_Shape_rec/bHLH_log10pVal_Pos.csv')
pVal.Pos[is.na(pVal.Pos)] <- 0
gplots::heatmap.2(as.matrix(as.data.frame(lapply(pVal.Pos, as.numeric))),dendrogram='none',
                  Rowv=FALSE, Colv=FALSE,trace='none',col = rev(heat.colors(12)), key.title = 'a')


####HT vs HT validation (use unfiltered mono_matrix)####
bHLH_motifs <- list()
for(i in 1:length(mono_motifs)){
  add <- list(name = mono_motifs[[i]]$name, matrix = mono_motifs[[i]]$monoMatrix)
  bHLH_motifs[[length(bHLH_motifs)+1]] <- add
}

posMotif <- 'P-1'

ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
PosAA <- c()
for(i in 1:length(bHLH_motifs)){
  PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
}
pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)

geneList <- strsplit(colnames(pos_matrix), '_')
geneNames <- c()
for(i in 1:length(geneList)){
  geneNames <- c(geneNames, geneList[[i]][1])
}
geneNames <- unique(geneNames)


P.1HT1 <- data.frame(placeHolder = c(0,0,0,0))
P.1HT2 <- data.frame(placeHolder = c(0,0,0,0))
for(i in 1:length(geneNames)){
  HTList <- grep(geneNames[i], colnames(pos_matrix))
  if(length(HTList) < 2){
    next()
  }
  P.1HT1 <- cbind.data.frame(P.1HT1, pos_matrix[,HTList[1]])
  P.1HT2 <- cbind.data.frame(P.1HT2, pos_matrix[,HTList[2]])
}

P.1HT1 <- P.1HT1[,-1]
P.1HT2 <- P.1HT2[,-1]

P.1HT1 <- log(P.1HT1)
P.1HT2 <- log(P.1HT2)
HT1 <- unlist(data.frame(apply(P.1HT1, 2, function(column) column - mean(column))))
HT2 <- unlist(data.frame(apply(P.1HT2, 2, function(column) column - mean(column))))

#HT1 <- unlist(data.frame( matrix2tetrahedron(P.1HT1)))
#HT2 <- unlist(data.frame(matrix2tetrahedron(P.1HT2)))

#HT1 <- log(unlist(data.frame(P.1HT1)))
#HT2 <- log(unlist(data.frame(P.1HT2)))

HTdata <- data.frame(HT1, HT2)

plot(HT1,HT2, pch = 19, xlab = 'Experiment 1 ¦¤¦¤G', ylab = 'Experiment 2 ¦¤¦¤G', main = 'Pair-wise experiment comparison', col = c('green','blue','orange','red'))
legend( x = "bottomright", legend = c('A', 'C', 'G', 'T'), col =c('green','blue','orange','red'), pch = 19, bty='n')

lm <- lm(HT2~HT1+0, HTdata)
summary(lm)
mean((HT2-HT1)^2)^(1/2)





####visualize tetrahedron####
AApos <- 13
posMotif <- 'P-1'

ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
PosAA <- c()
for(i in 1:length(bHLH_motifs)){
  PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
}
pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)
ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
PosAA <- c()
for(i in 1:length(bHLH_motifs)){
  PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
}
pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)
colnames(pos_matrix) <- PosAA
levels(as.factor(colnames(pos_matrix)))
global_plot <- plot_tetrahedron_resiColor(pos_matrix, size = 10,base_colors = c('green','blue','orange','red'), 
                                          resiColors = c('#cc9933', '#ff9933', '#cc99cc','#990000','#00ffff','#ffcc33'))
global_plot
htmlwidgets::saveWidget(as_widget(global_plot), "E:/Desktop/HB_Rotation_1/Base_Shape_rec/P-1_AA13_tetrahedron/tetrahedron.html")


####get predict eigen####
posMotif <- 'P-1'

ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
PosAA <- c()
for(i in 1:length(bHLH_motifs)){
  PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
}
pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)

geneList <- strsplit(colnames(pos_matrix), '_')
geneNames <- c()
for(i in 1:length(geneList)){
  geneNames <- c(geneNames, geneList[[i]][1])
}
geneNames <- unique(geneNames)

predList <- c()
TrueList <- c()
keyPos = c(13,14,5,8,26)
test.no <- length(geneNames)
toPs <- sample(1:length(geneNames), test.no)
for(a in 1:test.no){
  j <- toPs[a]
  pTo <- grep(geneNames[j], colnames(pos_matrix))[1]
  entryTo <- colnames(pos_matrix)[pTo]
  targetMatrix <- data.frame(pos_matrix[,entryTo])
  targetSequence <- bHLH_pbAlignment[bHLH_pbAlignment$name == entryTo,2]
  dddGList <- form.ddGFeatureList(bHLH_pbAlignment, bHLH_motifs, leaveOut = c(geneNames[j]), keyPos = keyPos)
  predictedMatrix <- predEigen(targetSequence, dddGList, motifPos = 'P-1', keyPos = keyPos, distType = 'Euclidean', singularity = 'one')
  predList <- c(predList, unlist(predictedMatrix))
  TrueList <- c(TrueList, unlist(data.frame(apply(log(targetMatrix), 2, function(column) column - mean(column)))))
}

dt <- data.frame(TrueList, predList)
plot(TrueList, predList, main = 'Pair-wise prediction test', xlab = 'Experiment ¦¤¦¤G', ylab = 'Predicted ¦¤¦¤G')
summary(lm(predList~TrueList, dt))
mean((TrueList[!is.na(predList)]-predList[!is.na(predList)])^2)^(1/2)
