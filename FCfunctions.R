write.bash <- function(script, outfile){
  cat("#!/bin/bash\n", file = outfile, sep = "\n", append = F)
  for(i in 1:length(script)){
    cat(script[i], file = outfile, sep = "\n", append = T)
  }
  return('bash script done, remember to dos2unix')
}

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

matrix2seq <- function(scoringMatrix){
  seq <- ''
  for(i in 1:ncol(scoringMatrix)){
    add <- rownames(scoringMatrix)[scoringMatrix[,i] == max(scoringMatrix[,i])]
    if(length(add) > 1){
      add <- 'N'
    }
    seq <- paste(seq, add, sep = '')
  }
  return(seq)
}

pymolOpenFiles <- function(RA.files){
  pymol <- c()
  for(i in 1:length(RA.files)){
    pymol[i] <- paste('load ', RA.files[i], sep = '')
  }
  pymol <- c(pymol,'hide all', 'select polymer.protein', 'show stick, sele', 'select polymer.nucleic', 'show line, sele', 'hide lines, hydrogen', 'hide stick, hydrogen')
  run.pymol(pml = pml, script = pymol)  
}


bpAA <- function(SNAPList){
  bp <- substr(SNAPList$bp.aa,1,2)
  AA <- substr(SNAPList$bp.aa,4,6)
  bpAA <- data.frame(bp, AA)
  return(bpAA)
}


extractPair.OriginAA <- function(
  SNAPList.AA,
  pml,
  SNAPdir,
  inspect.file,
  dssr = 'E:/x3dna-dssr.exe',
  cutoff = 4.5
){
  if(!file.exists(inspect.file)){
    dir.create(inspect.file)
  }
  Files <- c()
  for(i in 1:nrow(SNAPList.AA)){
    pairFile <- extractPair(testSNAP = SNAPList.AA[i,], pml = pml, SNAPFile = SNAPdir, 
                            pdbFile = inspect.file, tag = paste('_',SNAPList.AA$aa[i], sep = ''),
                            cutoff = cutoff)
    AAcode <- paste(substr(SNAPList.AA$aa[i],1,2), substr(SNAPList.AA$aa[i], 6, nchar(SNAPList.AA$aa[i])), sep = '')
    cmd <- paste(dssr, ' -i=', pairFile,' --frame-aa=', AAcode,' -o=',pairFile, sep = '')
    system(cmd)
    Files <- c(Files, pairFile)
  }
  return(Files)
}


screenSNAP.AA <- function(
  PDBdir,
  SNAPdir,
  pdbid ,
  hmmIndex,
  AApos,
  SNAP.update = TRUE,
  block = "pair/amino-acid interactions",
  cutoff = 4.5
){
  out <- data.frame(NULL)
  for(j in 1:length(pdbid)){
    snapName <- list.files(SNAPdir)[grep(pdbid[j], list.files(SNAPdir))]
    if(length(snapName) != 1 || SNAP.update){
      snap <- formSNAP(PDBdir, SNAPdir, pdbid[j], dssr, cutoff = cutoff)
    }else{
      snap <- paste(SNAPdir, '/', snapName, sep = '')
    }
    snap_df <- readSNAP(snap, block)
    indexPDB <- hmmIndex[hmmIndex$pdbid == pdbid[j],]
    for(i in 1:nrow(indexPDB)){
      chain <- indexPDB$chain[i]
      pos <- indexPDB[i,colnames(indexPDB) == paste('AA', AApos, sep = '')]
      add <- snap_df[substr(snap_df$aa,1,1)==chain,]
      add <- add[substr(add$aa,6,nchar(add$aa))==pos,]
      out <- rbind.data.frame(out,add)
    }
  }
  return(out)
}





readHmm <- function(hmmFile){
  hmmLines <- read.csv(hmmFile, header=FALSE, sep=";")
  hmmAlign <- read.csv(hmmFile,header = FALSE, quote="\"", comment = '#', sep = '')
  AlLen <- nchar(hmmAlign[1,2])
  #organize hmmLines
  GS <- substr(hmmLines[,1],1,4) == '#=GS'
  hmmLines <- hmmLines[GS,]
  
  hmmIndex <- data.frame(NULL)
  for(i in 1:length(hmmLines)){
    linetest <- hmmLines[i]
    split1 <- strsplit(linetest, '\\[')
    split2 <- strsplit(split1[[1]], '\\]')
    id <- split2[[1]][1]
    id <- strsplit(id, ' ')[[1]]
    id <- id[grep('/', id)]
    for(j in 4:length(split2)){
      addline <- c()
      entry <- split2[[j]][1]
      entry <- strsplit(entry, ', ')[[1]]
      pos <- strsplit(entry[3], '-')[[1]]
      addline <- c(addline, entry[1:2], pos, id)
      ind <- as.numeric(pos[1])
      seq <- hmmAlign[hmmAlign$V1==id,2]
      for(k in 1:AlLen){
        AA <- substr(seq, k,k)
        if(AA == '.' || AA == '-'){
          addline <- c(addline, NA)
        }else{
          addline <- c(addline, ind)
          ind <- ind+1
        }
      }
      hmmIndex <- rbind.data.frame(hmmIndex, addline)
    }
  }
  colnames(hmmIndex) <- c('pdbid', 'chain', 'start', 'end', 'alignment', paste('AA', c(1:AlLen), sep = ''))
  return(hmmIndex)
}


run.pymol <- function(pymol.dir = 'E:/pymol/pymol_app/pyMOLWin.exe',
                      pml = '~/pymolScript.pml',
                      script){
  cat('#run pymol \n', file = pml, append = F)
  for(i in 1:length(script)){
    cat(script[i],' \n', sep = '', file = pml, append = T)
  }
  cmd <- paste(pymol.dir, ' ', pml, sep = '')
  system(cmd, wait = TRUE)
  return('Pymol started')
}


AnchorOriginSNAP <- function(PDBdir,
                             SNAPdir,
                             anchorID,
                             line,
                             block = 'pair/amino-acid interactions',
                             pml,
                             inspect.file,
                             cutoff = 4.5
                             ){
  Anchor <- selectAnchor(PDBdir, SNAPdir, pdbid = anchorID, block = block, line = line, cutoff = cutoff)
  Adt <- rbind.data.frame(Anchor,Anchor)
  colnames(Adt) <- names(Anchor)
  Adt <- Adt[1,]
  ApdbFile <- extractPair(Adt, pml, SNAPFile = SNAPdir, pdbFile = inspect.file, block = block, cutoff = cutoff)
  Anchor.pdb <- bio3d::read.pdb(ApdbFile)
  Anchor.origin <- Anchor.pdb$atom[Anchor.pdb$atom$elety == 'CA',c('x','y','z')]
  Anchor.new <- c(Anchor, Anchor.origin)
  attr(Anchor.new, 'block') <- attr(Anchor,'block')
  attr(Anchor.new, 'cutoff') <- attr(Anchor, 'cutoff')
  Anchor <- Anchor.new
  return(Anchor)
}


ScreenOriginSNAP <- function(PDBdir, SNAPdir, pdbid, pml, Anchor, inspect.file, structure = F,
                             cutoff = 4.5,
                             threshold = c(5,50),
                             tag = '',
                             block = 'pair/amino-acid interactions',
                             pymol.dir = 'E:/pymol/pymol_app/pyMOLWin.exe',
                             wait.time = 3,
                             SNAP.update = T
                             ){
    if(SNAP.update){
      snap.file <- formSNAP(PDBdir = PDBdir, SNAPdir = SNAPdir, pdbid = pdbid, auxfile = T, cutoff = cutoff)
    }else{
      snap.file <- paste()
    }
    
    snap <- readSNAP(snap.file, block = attr(Anchor,'block'))
    pairs <- unique(snap$bp.aa)
    an.ch <- unlist(Anchor[c('Cx','Cy','Cz','Rx','Ry','Rz')])
    if(length(unique(pairs)) < 1){
      return()
    }
    snap.xyz <- data.frame(NULL)
    for(i in 1:length(unique(pairs))){
      snap.pair <- snap[snap$bp.aa == pairs[i],]
      for(j in 1:nrow(snap.pair)){
        sa.ch <- snap.pair[j,c('Cx','Cy','Cz','Rx','Ry','Rz')]
        distOR <- CAdist(as.numeric(an.ch),as.numeric(sa.ch))
        if(distOR[1] <= threshold[1] && distOR[2] <= threshold[2]){
          add <- c(unlist(snap.pair[j,]),distO = distOR[1], distR = distOR[2])
          snap.xyz <- rbind.data.frame(snap.xyz, add)
          pdb.file <- paste(PDBdir, '/',pairs[i], '.pdb', sep = '')
          colnames(snap.xyz) <- names(add)
          if(structure){
            cat('#save state \n', file = pml, append = F)
            cat('load ', pdb.file,' \n', sep = '', file = pml, append = T)
            outfile <- paste(inspect.file, '/', snap.pair$id[1], '_', pairs[i], '_', j,'_', tag,'_snap', '.pdb', sep = '')
            cat('save ', outfile,', state=', j, ' \n', sep = '', file = pml, append = T)
            cat('quit \n', file = pml, append = T)
            cmd <- paste(pymol.dir, ' -qc ', pml, sep = '')
            system(cmd, wait = TRUE)
            Sys.sleep(wait.time)
          }
        }
      }
    }
    options(warn=-1)
    files <- list.files(PDBdir)
    deleFiles <- files[!complete.cases(as.numeric(substr(files,1,1)))]
    file.remove(paste(PDBdir, '/', deleFiles, sep = ''))
    options(warn=0)
  return(snap.xyz)
}


CAdist <- function(anchor, sample){
  distO <- dist(rbind(anchor[1:3], sample[1:3]))
  distR <- (rotationDiff(anchor[4], sample[4])^2 + rotationDiff(anchor[5], sample[5])^2 + 
              rotationDiff(anchor[6], sample[6])^2)^(1/2)
  return(c(distO,distR))
}

rotationDiff <- function(R1, R2){
  P1 <- R1-R2
  if(P1 > 180){
    P1 <- -(360-P1)
  }else if(P1 < -180){
    P1 <- 360+P1
  }
  return(P1)
}



TR2matrix <- function(TR, R.first = T){ #colnames have to be Tx, Ty, Tz, Rx, Ry, Rz
  out <- matrix(nrow = 4, ncol = 3, data = c(0,0,0,1,0,0,0,1,0,0,0,1), byrow = T)
  colnames(out) <- c('x', 'y', 'z')
  rownames(out) <- c('origin', 'Xaxis', 'Yaxis', 'Zaxis')
  cols <- names(TR)
  xTh <- as.numeric(TR[4])/180*pi
  yTh <- as.numeric(TR[5])/180*pi
  zTh <- as.numeric(TR[6])/180*pi
  Tx <- as.numeric(TR[1])
  Ty <- as.numeric(TR[2])
  Tz <- as.numeric(TR[3])
  RxM <- matrix(nrow = 3, ncol = 3, data = c(1,0,0,0,cos(xTh), -sin(xTh), 0, sin(xTh), cos(xTh)), byrow = T)
  RyM <- matrix(nrow = 3, ncol = 3, data = c(cos(yTh),0,sin(yTh),0,1, 0, -sin(yTh), 0, cos(yTh)), byrow = T)
  RzM <- matrix(nrow = 3, ncol = 3, data = c(cos(zTh), -sin(zTh), 0,sin(zTh), cos(zTh), 0,0,0,1 ), byrow = T)
  if(R.first){
    for(i in 2:nrow(out)){
      out[i,] <- RzM%*%RyM%*%RxM%*%out[i,]
    }
    out[1,] <- out[1,] + Tx*out[2,] + Ty*out[3,] + Tz*out[4,]
  }else{
    out[1,] <- out[1,] + Tx*out[2,] + Ty*out[3,] + Tz*out[4,]
    for(i in 2:nrow(out)){
      out[i,] <- RzM%*%RyM%*%RxM%*%out[i,]
    }
  }
  return(out)
}

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

TRmatrixDist <- function(TRM1, TRM2){
  OriginDist <- dist(rbind.data.frame(TRM1[1,], TRM2[1,]))
  A1 <- angle(TRM1[2,], TRM2[2,])
  A2 <- angle(TRM1[3,], TRM2[3,])
  A3 <- angle(TRM1[4,], TRM2[4,])
  AngleDist <- (A1^2 + A2^2 + A3^2)^(1/2)
  return(c(OriginDist,AngleDist))
}

extractPair <- function(testSNAP,
                        pml,
                        SNAPFile,
                        pdbFile,
                        block = 'pair/amino-acid interactions',
                        pymol.dir = 'E:/pymol/pymol_app/pyMOLWin.exe',
                        wait.time = 3,
                        cutoff = 4.5,
                        tag = ''){
  out <- c()
  for(i in 1:nrow(testSNAP)){
    snap <- formSNAP(PDBdir = pdbFile, 
                     SNAPdir = SNAPFile, 
                     pdbid = testSNAP$id[i], auxfile = T, cutoff = cutoff)
    snap_df <- as.data.frame(readSNAP(snap,block))
    pair <- testSNAP$bp.aa[i]
    snap_df <- snap_df[snap_df$bp.aa == pair,]
    state <- intersect(which(as.numeric(snap_df$Tdst) == as.numeric(testSNAP$Tdst[i])), which(as.numeric(snap_df$Rdst) == as.numeric(testSNAP$Rdst[i])))
    pairFile <- paste(pair, '.pdb', sep = '')
    cat('#save state \n', file = pml, append = F)
    cat('load ', pdbFile, '/',pairFile,' \n', sep = '', file = pml, append = T)
    outfile <- paste(pdbFile, '/', testSNAP$id[i], '_', pair, '_', state, tag, '.pdb', sep = '')
    cat('save ', pdbFile, '/', testSNAP$id[i], '_', pair, '_', state, tag, '.pdb, state=', state, ' \n', sep = '', file = pml, append = T)
    cat('quit \n', file = pml, append = T)
    cmd <- paste(pymol.dir, ' -qc ', pml, sep = '')
    system(cmd, wait = TRUE)
    Sys.sleep(wait.time)
    snap <- formSNAP(PDBdir = pdbFile,
                     SNAPdir = SNAPFile, 
                     pdbid = testSNAP$id[i], auxfile = F, cutoff = cutoff)
    out <- c(out, outfile)
  }
  options(warn=-1)
  files <- list.files(pdbFile)
  deleFiles <- files[!complete.cases(as.numeric(substr(files,1,1)))]
  file.remove(paste(pdbFile, '/', deleFiles, sep = ''))
  options(warn=0)
  return(out)
}



addBaseLink <- function(testScreen, PDBdir){
  baseLink <- c()
  for(i in 1:nrow(testScreen)){
    x <- tryCatch({
      a1 <- readPosition(bio3d::read.pdb(download.pdb(testScreen$pdb[i], PDBdir,F)), baseChain = testScreen$baseChain[i], baseResi = testScreen$baseNO[i], 
                       AAChain = testScreen$AAChain[i], AAResi = testScreen$AANo[i],
                       baseAnchor = c('N1'), AAAnchor = c('CA'))
    
      a2 <- readPosition(bio3d::read.pdb(download.pdb(testScreen$pdb[i], PDBdir,F)), baseChain = testScreen$base2Chain[i], baseResi = testScreen$base2NO[i], 
                         AAChain = testScreen$AAChain[i], AAResi = testScreen$AANo[i],
                         baseAnchor = c('N1'), AAAnchor = c('CA'))
    }, error = function(e) e) 
    if(length(x) == 2){
      baseLink <- c(baseLink,-1)
    }else if(a1 < a2){
      baseLink <- c(baseLink, testScreen$baseType[i])
    }else{
      baseLink <- c(baseLink, testScreen$base2Type[i])
    }
  }
  return(cbind.data.frame(testScreen, baseLink))
}

getDNAPos <- function(string){
  posList <- c()
  posList <- c(posList, gregexpr('A', string)[[1]])
  posList <- c(posList, gregexpr('T', string)[[1]])
  posList <- c(posList, gregexpr('C', string)[[1]])
  posList <- c(posList, gregexpr('G', string)[[1]])
  if(length(posList) > 0){
    return(max(posList))
  }else{
    return(-1)
  }
  
}

SNAP2screenList <- function(SNAPList){
  out <- data.frame(NULL)
  if(length(grep('nt2', colnames(SNAPList))) > 0){
    for(i in 1:nrow(SNAPList)){
      pdb <- SNAPList$id[i]
      base1 <- strsplit(SNAPList$nt1[i],'\\.')[[1]]
      baseChain <- base1[1]
      pos <- getDNAPos(base1[2])
      baseType <- substr(base1[2], 1, pos)
      baseNO <- substr(base1[2], pos+1, nchar(base1[2]))
      base1 <- strsplit(SNAPList$nt2[i],'\\.')[[1]]
      base2Chain <- base1[1]
      pos <- getDNAPos(base1[2])
      base2Type <- substr(base1[2], 1, pos)
      base2NO <- substr(base1[2], pos+1, nchar(base1[2]))
      base1 <- strsplit(SNAPList$aa[i],'\\.')[[1]]
      AAChain <- base1[1]
      AAType <- gsub('[[:digit:]]', '', base1[2])
      AANo <- gsub(AAType, '', base1[2])
      
      add <- data.frame(pdb, baseChain, baseNO, baseType, base2Chain, base2NO, base2Type,
                        AAChain, AANo, AAType, maxDiff = as.numeric(SNAPList$dist[i]))
      out <- rbind.data.frame(out, add)
    }
  }else{
    for(i in 1:nrow(SNAPList)){
      pdb <- SNAPList$id[i]
      base1 <- strsplit(SNAPList$nt[i],'\\.')[[1]]
      baseChain <- base1[1]
      baseType <- gsub('[[:digit:]]', '', base1[2])
      baseNO <- gsub(baseType, '', base1[2])
      base1 <- strsplit(SNAPList$aa[i],'\\.')[[1]]
      AAChain <- base1[1]
      AAType <- gsub('[[:digit:]]', '', base1[2])
      AANo <- gsub(AAType, '', base1[2])
      
      add <- data.frame(pdb, baseChain, baseNO, baseType, 
                        AAChain, AANo, AAType, maxDiff = as.numeric(SNAPList$dist[i]))
      out <- rbind.data.frame(out, add)
    }
  }
  return(out)
}



pymolInspect <- function(screenList, projectFile, outFile, scriptFile, outTag = '', base.no = 1, adj = F){
  cat('#screenList for pyMol \n', file = scriptFile, append = F)
  for(i in 1:nrow(screenList)){
    pdbname <- screenList$pdb[i]
    cat('load ', projectFile, '/',pdbname,'.pdb \n', sep = '', file = scriptFile, append = T)
    if(base.no == 1){
      cat('select (chain ',screenList$baseChain[i], ' and resi ', screenList$baseNO[i], ') or (chain ',
          screenList$AAChain[i], ' and resi ', screenList$AANo[i], ') \n', sep = '', file = scriptFile, append = T)
    }else if(base.no == 2){
      if(adj){
        cat('select (chain ',screenList$baseChain[i], ' and resi ', as.numeric(screenList$baseNO[i]), '-', 
            as.numeric(screenList$baseNO[i])+1, ') or (chain ',
            screenList$AAChain[i], ' and resi ', screenList$AANo[i], ') or (chain ', screenList$base2Chain[i], ' and resi ',
            as.numeric(screenList$base2NO[i]), '-', as.numeric(screenList$base2NO[i])+1,
             ') \n', sep = '', file = scriptFile, append = T)
      }else{
        cat('select (chain ',screenList$baseChain[i], ' and resi ', screenList$baseNO[i], ') or (chain ',
            screenList$AAChain[i], ' and resi ', screenList$AANo[i], ') or (chain ', screenList$base2Chain[i], ' and resi ', 
            screenList$base2NO[i], ') \n', sep = '', file = scriptFile, append = T)
      }
    }else{
      print('error in base.no')
      return()
    }
    
    cat('create obj01, sele \n', sep = '', file = scriptFile, append = T)
    cat('delete ', pdbname, ' \n', sep = '', file = scriptFile, append = T)
    cat('show sticks, obj01 \n', file = scriptFile, append = T)
    cat('save ', outFile, '/', pdbname, '_', outTag, '.pdb \n', sep = '', file = scriptFile, append = T)
    cat('delete all \n', file = scriptFile, append = T)
  }
  cat('quit \n', file = scriptFile, append = T)
}



screenPosition <- function(pdbList,
                           projectFile,
                           threshold = 0.5,
                           anchorMatrix){
  out <- data.frame(NULL)
  for(i in 1:length(pdbList)){
    tryCatch({
      pdbname <- pdbList[i]
      pdb <- download.pdb(pdbname, projectFile = projectFile, as.file = FALSE)
      pdb <- read.pdb(pdb)
      atom <- pdb$atom
      baseAtoms <- data.frame(NULL)
      AAAtoms <- data.frame(NULL)
      for(j in 1:nrow(atom)){
        Line <- atom[j,]
        if(Line$resid == 'DA' || Line$resid == 'DC' || Line$resid == 'DG'|| Line$resid == 'DT'){
          baseAtoms <- rbind.data.frame(baseAtoms, Line)
        }else if(Line$resid == 'HOH'){
          next()
        }else{
          AAAtoms <- rbind.data.frame(AAAtoms, Line)  
        }
      }
      baseUse <- data.frame(NULL)
      for(k in 1:nrow(anchorMatrix)){
        baseAdd <- baseAtoms[grep(rownames(anchorMatrix)[k], baseAtoms$elety),]
        baseAdd <- baseAdd[nchar(baseAdd$elety) == nchar(rownames(anchorMatrix)[k]),]
        baseUse <- rbind.data.frame(baseUse,baseAdd)
      }
      basePairs <- unique(baseUse[,c('chain','resno')])
      AAUse <- data.frame(NULL)
      for(k in 1:ncol(anchorMatrix)){
        AAAdd <- AAAtoms[grep(colnames(anchorMatrix)[k], AAAtoms$elety),]
        AAAdd <- AAAdd[nchar(AAAdd$elety) == nchar(colnames(anchorMatrix)[k]),]
        AAUse <- rbind.data.frame(AAUse,AAAdd)
      }
      AAPairs <- unique(AAUse[,c('chain','resno')])
      for(k in 1:nrow(basePairs)){
        baseNow <- baseUse[baseUse[,'chain'] == basePairs[k,'chain'],]
        baseNow <- baseNow[baseNow[,'resno'] == basePairs[k,'resno'],]
        if(nrow(baseNow) < nrow(anchorMatrix)){
          next()
        }
        for(l in 1:nrow(AAPairs)){
          AANow <- AAUse[AAUse[,'chain'] == AAPairs[l,'chain'],]
          AANow <- AANow[AANow[,'resno'] == AAPairs[l,'resno'],]
          if(nrow(AANow) < ncol(anchorMatrix)){
            next()
          }
          fit <- 0
          maxDiff <- -Inf
          for(A1 in 1:nrow(anchorMatrix)){
            for(A2 in 1:ncol(anchorMatrix)){
              true <- anchorMatrix[A1,A2]
              XYZbase <- baseNow[grep(rownames(anchorMatrix)[A1],baseNow$elety),c('x','y','z')]
              XYZAA <- AANow[grep(colnames(anchorMatrix)[A2],AANow$elety),c('x','y','z')]
              if(nrow(XYZbase) != 1 || nrow(XYZAA) != 1){
                next
              }
              test <- dist(rbind(XYZbase,XYZAA))
              if(abs(true - test) < threshold){
                fit <- fit + 1
                if(abs(true-test) > maxDiff){
                  maxDiff <- abs(true-test)
                }
              }
            }
          }
          if(fit == nrow(anchorMatrix)*ncol(anchorMatrix)){
            add <- data.frame(pdb = pdbList[i], baseChain = baseNow$chain[1], baseNO = baseNow$resno[1], baseType = baseNow$resid[1],
                              AAChain = AANow$chain[1], AANo = AANow$resno[1], AAType = AANow$resid[1], maxDiff = as.numeric(maxDiff))
            out <- rbind.data.frame(out, add)
            print('Base-AA pair found')
          }
        }
      }
      #print(paste('screened over pdb: ', pdbname, '; pdb number:', i, sep = ''))
    }, error = function(e) e)
  }
  out <- unique(out)
  return(out)
}

#function




readPosition <- function(pdb, 
                         baseChain = 'B',
                         baseResi = 9,
                         AAChain = 'A',
                         AAResi = 217,
                         baseAnchor = c('P','N9','C1','O4'),
                         AAAnchor = c('CA','CZ')){
  atom <- pdb$atom
  base <- atom[atom$chain == baseChain,]
  base <- base[base$resno == baseResi,]
  AA <- atom[atom$chain == AAChain,]
  AA <- AA[AA$resno == AAResi,]
  outMatrix <- matrix(nrow = length(baseAnchor), ncol = length(AAAnchor), data = NA)
  rownames(outMatrix) <- baseAnchor
  colnames(outMatrix) <- AAAnchor
  for(i in 1:length(baseAnchor)){
    baseAtom <- base[grep(baseAnchor[i], base$elety),]
    baseAtom <- baseAtom[nchar(baseAtom$elety) == min(nchar(baseAtom$elety)),]
    baseXYZ <- c(baseAtom$x,baseAtom$y,baseAtom$z)
    for(j in 1:length(AAAnchor)){
      AAAtom <- AA[grep(AAAnchor[j], AA$elety),]
      AAAtom <- AAAtom[nchar(AAAtom$elety) == min(nchar(AAAtom$elety)),]
      AAXYZ <- c(AAAtom$x,AAAtom$y,AAAtom$z)
      addDist <- dist(rbind(baseXYZ,AAXYZ))
      outMatrix[i,j] <- addDist
    }
  }
  return(outMatrix)
}


formSNAP <- function(PDBdir, SNAPdir = PDBdir, pdbid, dssr = 'E:/x3dna-dssr.exe', auxfile = F, cutoff = 4.5){
  wd <- getwd()
  setwd(PDBdir)
  proteinStructureBoost::download.pdb(pdbid, projectFile = PDBdir,as.file = FALSE)
  pdbfile <- paste(pdbid,'.pdb',sep = '')
  SNAPfile <- paste(SNAPdir, '/', pdbid, '_SNAP.txt', sep = '')
  if(auxfile){
    cmd <- paste(dssr, ' snap -i=', pdbfile, ' --c-alpha --pair-order=as-is --cutoff=',cutoff,' --auxfile > ', SNAPfile, sep = '')
  }else{
    cmd <- paste(dssr, ' snap -i=', pdbfile, ' --c-alpha --pair-order=as-is --cutoff=',cutoff,' > ', SNAPfile, sep = '')
  }
  shell(cmd)
  setwd(wd)
  return(SNAPfile)
}

readSNAP <- function(SNAPfile, block = 'nucleotide/amino-acid interaction'){
  options(warn=-1)
  data <- read.table(file = SNAPfile,sep='\n',fill=TRUE,header = F)
  start <- grep(block, data$V1) + 2
  if(length(start) != 1){
    print('error in block value')
    return()
  }
  end <- min(grep('\\*', data$V1)[grep('\\*', data$V1) > start]) - 1
  if(length(end) == 0){
    end <- nrow(data)
  }
  trim <- strsplit(gsub('\\s+', '\t', data[start:end,]), split = '\t')
  if(length(trim) > 999){
    for(tt in 1000:length(trim)){
      trim[[tt]] <- c('',trim[[tt]])
    }
  }
  trim <- as.data.frame(trim) 
  trim <- t(trim)
  trim <- trim[,3:ncol(trim)]
  colnames(trim) <- strsplit(gsub('\\s+', '\t', data[start-1,]), split = '\t')[[1]][-1]
  colnames(trim)[2] <- 'bp.aa'
  rownames(trim) <- 1:nrow(trim) 
  options(warn=0)
  return(as.data.frame(trim))
}

selectAnchor <- function(PDBdir, SNAPdir = PDBdir, pdbid, dssr = 'E:/x3dna-dssr.exe', line = NA, block = 'nucleotide/amino-acid interaction', cutoff = 4.5){
  snap <- formSNAP(PDBdir, SNAPdir, pdbid, dssr, cutoff = cutoff)
  snap_df <- readSNAP(snap, block)
  if(is.na(line)){
    print(snap_df)
    anchorLine <- my.name <- readline(prompt="Enter anchor interaction line number: ")
  }else{
    anchorLine <- line
  }
  anchorLine <- snap_df[anchorLine,]
  attr(anchorLine, 'block') <- block
  attr(anchorLine, 'cutoff') <- cutoff
  return(anchorLine)
}

screenSNAP <- function(PDBdir, SNAPdir = PDBdir, anchorLine, pdbid,  dssr = 'E:/x3dna-dssr.exe', 
                       threshold = c(5,1), method = 'matrix'){ #avaliable methods: raw (8 thresh), abs (8 thresh), matrix (2 thresh)
  block <- attr(anchorLine, 'block')
  cutoff <- attr(anchorLine, 'cutoff')
  out <- data.frame(NULL)
  snapName <- list.files(SNAPdir)[grep(pdbid, list.files(SNAPdir))]
  if(length(snapName) != 1){
    snap <- formSNAP(PDBdir, SNAPdir, pdbid, dssr, cutoff = cutoff)
  }else{
    snap <- paste(SNAPdir, '/', snapName, sep = '')
  }
  snap_df <- readSNAP(snap, block)
  if(is.null(snap_df)){
    attr(out, 'is.dep') <- 1
    return(out)
  }
  anchorCoord <- anchorLine[c('Tdst', 'Rdst', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz')]
  snapCoords <- snap_df[,c('Tdst', 'Rdst', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz')]
  for(i in 1:nrow(snap_df)){
    if(method == 'raw'){
      dist <- SNAP8dDist.abs(anchorCoord, snapCoord = snapCoords[i,], threshold = threshold)
    }else if(method == 'abs'){
      dist <- SNAP8dDist(anchorCoord, snapCoord = snapCoords[i,], threshold = threshold)
    }else if(method == 'matrix'){
      dist <- SNAPMatrixDist(anchorCoord, snapCoord = snapCoords[i,], threshold = threshold)
      if(as.numeric(anchorCoord[2])*as.numeric(snapCoords[i,2]) < 0){
        dist <- Inf
      }
    }
    if(dist <= 1){
      dist <- as.numeric(dist)
      line <- c(snap_df[i,], dist = dist)
      line <- t(as.data.frame(line))
      out <- rbind.data.frame(out, line)
    }
  }
  return(out)
}

SNAP8dDist <- function(anchorCoord, snapCoord, threshold = c(Inf, Inf, 1,1,1,Inf,Inf,Inf)){
  return(sum(((as.numeric(anchorCoord) - as.numeric(snapCoord))/threshold)^2)^(1/2))
}

SNAP8dDist.abs <- function(anchorCoord, snapCoord, threshold = c(Inf, Inf, 1,1,1,Inf,Inf,Inf)){
  return(sum(((abs(as.numeric(anchorCoord)) - abs(as.numeric(snapCoord)))/threshold)^2)^(1/2))
}

SNAPMatrixDist <- function(anchorCoord, snapCoord, threshold = c(5,0.5)){
  TRM1 <- TR2matrix(anchorCoord[3:8])
  TRM2 <- TR2matrix(snapCoord[3:8])
  dist <- TRmatrixDist(TRM1, TRM2)
  if(dist[1] <= threshold[1] && dist[2] <= threshold[2]){
    return(dist[1]/threshold[1])
  }else{
    return(Inf)
  }
}

SNAPsummary <- function(data){
  SNAP.summary <- matrix(nrow = 3, ncol = 8, data = 0)
  rownames(SNAP.summary) <- c('mean', 'median', 'sd')
  colnames(SNAP.summary) <- c('Tdst', 'Rdst', 'Tx', 'Ty', 'Tz', 'Rx', 'Ry', 'Rz')
  for(i in 1:length(colnames(SNAP.summary))){
    SNAP.summary[1,i] <- mean(as.numeric(data[,colnames(SNAP.summary)[i]]))
    SNAP.summary[2,i] <- median(as.numeric(data[,colnames(SNAP.summary)[i]]))
    SNAP.summary[3,i] <- sd(as.numeric(data[,colnames(SNAP.summary)[i]]))
  }
  return(SNAP.summary)
}



JSON2Matrix <- function(JSON_text, expo = exp(1)){
  library(RJSONIO)
  scoring <- fromJSON(JSON_text)
  scoring_matrix <- scoring$coefficients$bindingModes[[1]]$mononucleotide
  DN_seq <- strsplit(scoring$modelSettings$letterOrder, "")[[1]]
  N_pos <- scoring$modelSettings$bindingModes[[1]]$size
  if(expo < 0){
    out_matrix <- matrix(data = scoring_matrix, nrow = 4, ncol = N_pos, byrow = FALSE)
  }else{
    out_matrix <- matrix(data = expo^scoring_matrix, nrow = 4, ncol = N_pos, byrow = FALSE)
  }
  colnames(out_matrix) <- paste('P',c(1:N_pos),sep = '')
  rownames(out_matrix) <- DN_seq
  return(out_matrix)
}


normalize_sum1 <- function(mat, rows = TRUE){
  if(rows){
    mat <- t(mat)
  }
  for(i in 1:ncol(mat)){
    sum <- sum(mat[,i])
    for(j in 1:nrow(mat)){
      mat[j,i] <- mat[j,i] / sum
    }
  }
  if(rows){
    mat <- t(mat)
  }
  return(mat)
}


plot_tetrahedron <- function(JSON_matrix, base_colors = c('green','blue','orange','red'), size = 5, axis = FALSE){
  library(plotly)
  #transform to 3d space
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  df <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  AAtype <- rownames(df)
  df <- as.data.frame(df)
  df <- cbind.data.frame(df, name = AAtype)
  #create tetrahedron
  tetrahedron <-  matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1), nrow = 8, ncol = 3, byrow = TRUE)
  tetrahedron <- as.data.frame(tetrahedron)
  bases <- matrix(data = c(1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  bases <- as.data.frame(bases)
  #set parameters
  base_names <- rownames(JSON_matrix)
  len <- nrow(df)
  #set axis
  axx <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  #plot 3d graph
  plot3D <- plot_ly()%>%
    add_trace(tetrahedron, x = tetrahedron[,1],
              y=tetrahedron[,2],z=tetrahedron[,3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I('black'),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
              type = 'scatter3d',
              mode = 'text',
              text = df[1:len, 4],
              textposition = 'top',
              name = 'label',
              textfont = list(color = 'grey', size = size),
              opacity = 0.5)
  if(!axis){
    plot3D <- plot3D%>%layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx))
  }
  
  return(plot3D)
}

plot_tetrahedron_resiColor <- function(JSON_matrix, base_colors = c('green','blue','orange','red'), size = 5, axis = FALSE, resiColors = c('black', 'pink','magenta')){
  library(plotly)
  #transform to 3d space
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  df <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  AAtype <- rownames(df)
  df <- as.data.frame(df)
  df <- cbind.data.frame(df, name = AAtype)
  colLen <- length(resiColors)
  colorNum <- as.factor(df$name)
  resiCol <- c()
  for(i in 1:length(colorNum)){
    index <- as.numeric(colorNum[i])%%colLen
    if(index == 0){
      index <- colLen
    }
    resiCol <- c(resiCol, resiColors[index])
  }
  
  #create tetrahedron
  tetrahedron <-  matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1), nrow = 8, ncol = 3, byrow = TRUE)
  tetrahedron <- as.data.frame(tetrahedron)
  bases <- matrix(data = c(1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  bases <- as.data.frame(bases)
  #set parameters
  base_names <- rownames(JSON_matrix)
  len <- nrow(df)
  #set axis
  axx <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  #plot 3d graph
  plot3D <- plot_ly()%>%
    add_trace(tetrahedron, x = tetrahedron[,1],
              y=tetrahedron[,2],z=tetrahedron[,3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size/2, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size/2, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size/2, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size/2, opacity = 1),
              name = base_names[4])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I(resiCol),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries',
              showlegend = F)
  
  
  if(!axis){
    plot3D <- plot3D%>%layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx))
  }
  
  return(plot3D)
}


find_binding_site <- function(JSON_matrix,rec_seq,threshold = 0.5){
  #set output data structure
  out <- list(start_pos = c(), high_scores = c(), scores = c(), max_pos = 0)
  #normalize JSON_matrix
  norm_matrix <- normalize_sum1(JSON_matrix, rows = FALSE)
  len_motif <- ncol(norm_matrix)
  len_rec_seq <- nchar(rec_seq)
  for(i in 1:(len_motif - len_rec_seq + 1)){
    score <- 0
    for(j in 1:len_rec_seq){
      if(substr(rec_seq,j,j) == 'N'){
        score <- score + threshold
      }else{
        score <- score + norm_matrix[substr(rec_seq,j,j),(i+j-1)]
      }
    }
    score <- score/len_rec_seq
    if(score >= threshold){
      out$start_pos <- c(out$start_pos, i)
      out$high_scores <- c(out$high_scores, score)
    }
    out$scores <- c(out$scores, score)
    out$max_pos <- grep(max(out$scores), out$scores)
  }
  return(out)
}

#input a list of list objects with name = gene name & matrix = JSON matrix with unified colnames & rownames. pos enters as colnames
gene2pos <- function(motifs, pos = 'P1', nrow = 4){
  ngene <- length(motifs)
  out <- matrix(nrow = nrow, ncol = ngene, data = 0)
  rownames(out) <- rownames(motifs[[1]]$matrix)
  names <- c()
  for(i in 1:ngene){
    motif <- motifs[[i]]
    names <- c(names, motif$name)
    out[,i] <- motif$matrix[,pos]
  }
  colnames(out) <- names
  return(out)
}

plot_4Graph <- function(JSON_matrix, base_colors = c('green','blue','orange','red'), size = 5){
  library(plotly)
  #transform to 3d space
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  df <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  df <- as.data.frame(df)
  df <- cbind.data.frame(df, name = rownames(df))
  #create tetrahedron
  tetrahedron <-  matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1,1,1,1,-1,1,-1,1,-1,-1,-1,-1,1), nrow = 8, ncol = 3, byrow = TRUE)
  tetrahedron <- as.data.frame(tetrahedron)
  bases <- matrix(data = c(1,1,1, 1,-1,-1, -1,1,-1, -1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  bases <- as.data.frame(bases)
  #set parameters
  base_names <- rownames(JSON_matrix)
  len <- nrow(df)
  #set axis
  axx <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  #set cemera
  
  #plot 3d graph A
  plotA <- plot_ly(scene = 'scene1')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(2,3,4,2),], x = tetrahedron[c(2,3,4,2),1],
              y=tetrahedron[c(2,3,4,2),2],z=tetrahedron[c(2,3,4,2),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I('black'),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
              type = 'scatter3d',
              mode = 'text',
              text = df[1:len, 4],
              textposition = 'top',
              name = 'label',
              textfont = list(color = 'grey', size = size*2),
              opacity = 0.5)
  
  #plot 3d graph C
  plotC <- plot_ly(scene = 'scene2')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(1,3,4,1),], x = tetrahedron[c(1,3,4,1),1],
              y=tetrahedron[c(1,3,4,1),2],z=tetrahedron[c(1,3,4,1),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I('black'),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(bases, x = as.numeric(bases[4,1]), y = as.numeric(bases[4,2]), z = as.numeric(bases[4,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[4])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
              type = 'scatter3d',
              mode = 'text',
              text = df[1:len, 4],
              textposition = 'top',
              name = 'label',
              textfont = list(color = 'grey', size = size*2),
              opacity = 0.5)
  
  #plot 3d graph T
  plotT <- plot_ly(scene = 'scene4')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(1,2,3,1),], x = tetrahedron[c(1,2,3,1),1],
              y=tetrahedron[c(1,2,3,1),2],z=tetrahedron[c(1,2,3,1),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I('black'),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[3]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
              type = 'scatter3d',
              mode = 'text',
              text = df[1:len, 4],
              textposition = 'top',
              name = 'label',
              textfont = list(color = 'grey', size = size*2),
              opacity = 0.5)
  
  #reverse G plot
  tetra_trans_matrix <- matrix(data = c(1,-1,-1,1,1,1,-1,-1,1,-1,1,-1), nrow = 4, ncol = 3, byrow = TRUE)
  df2 <- normalize_sum1(t(JSON_matrix)) %*% tetra_trans_matrix
  df2 <- as.data.frame(df2)
  df <- cbind.data.frame(df2, name = rownames(df2))
  
  plotG2 <- plot_ly(scene = 'scene3')%>%
    #layout(scene = list(xaxis = axx, yaxis = axx, zaxis = axx, camera = eyepos))%>%
    add_trace(tetrahedron[c(1,2,3,1),], x = tetrahedron[c(1,2,3,1),1],
              y=tetrahedron[c(1,2,3,1),2],z=tetrahedron[c(1,2,3,1),3],color = I('blue'),
              type = 'scatter3d',
              mode = 'lines',
              line = list(width = 1),
              opacity = 0.5,
              name = 'tetrahedron')%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3], color = I('black'),
              type = 'scatter3d',
              mode = 'markers',
              size = size,
              opacity = 0.8,
              name = 'Entries')%>%
    add_trace(bases, x = as.numeric(bases[1,1]), y = as.numeric(bases[1,2]), z = as.numeric(bases[1,3]), color = I(base_colors[2]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[1])%>%
    add_trace(bases, x = as.numeric(bases[2,1]), y = as.numeric(bases[2,2]), z = as.numeric(bases[2,3]), color = I(base_colors[1]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[2])%>%
    add_trace(bases, x = as.numeric(bases[3,1]), y = as.numeric(bases[3,2]), z = as.numeric(bases[3,3]), color = I(base_colors[4]),
              type = 'scatter3d',
              mode = 'markers',
              marker = list(size = size, opacity = 1),
              name = base_names[3])%>%
    add_trace(df, x = df[1:len,1], y=df[1:len,2],z=df[1:len,3],
              type = 'scatter3d',
              mode = 'text',
              text = df[1:len, 4],
              textposition = 'top',
              name = 'label',
              textfont = list(color = 'grey', size = size*2),
              opacity = 0.5)
  
  plot4 <- subplot(plotA, plotC, plotG2, plotT)%>%
    layout(scene = list(domain=list(x=c(0,0.5),y=c(0,0.5)),
                        xaxis=axx, yaxis=axx, zaxis=axx,
                        camera = list(eye = list(x = 1, y = 1, z = 1)),
                        aspectmode='cube'),
           scene2 = list(domain=list(x=c(0.25,0.75),y=c(0.12,0.62)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         camera = list(eye = list(x = 1, y = -1, z = -1)),
                         aspectmode='cube'),
           scene3 = list(domain=list(x=c(0.25,0.75),y=c(0.43,0.93)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         camera = list(eye = list(x = -1, y = -1, z = 1)),
                         aspectmode='cube'),
           scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)),
                         xaxis=axx, yaxis=axx, zaxis=axx,
                         camera = list(eye = list(x = -1, y = -1, z = 1)),
                         aspectmode='cube'),
           showlegend = FALSE)
  return(plot4)
}


matrix2tetrahedron <- function(matrix){
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  return(normalize_sum1(t(matrix)) %*% tetra_trans_matrix)
}

mononucleotide_logo <- function(pwm, type="energy", axes=TRUE, reverse=FALSE,
                                labels=TRUE) {

  pwm <- apply(pwm, 2, function(column) column - mean(column))

  
  if (type == "info") {
    pwm <- exp(pwm)
  }
  
  if (type == "prob") {
    pwm <- apply(exp(pwm), 2, function(column) column / sum(column))
  }
  
  if (reverse) {
    pwm <- pwm[, rev(seq_len(ncol(pwm)))]
    
    if (is.null(dim(pwm))) {
      pwm <- matrix(pwm, length(pwm), 1)
    }
    
    rownames(pwm) <- letter_complement
  }
  
  if (type == "energy") {
    plot <- ggseqlogo(pwm, method = "c", font = font, col_scheme = col_scheme) +
      list(labs(x = NULL, y = NULL),
           annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                    alpha = 0.5, fill = "white"),
           scale_y_continuous(breaks = pretty_breaks()),
           geom_hline(yintercept = 0), alltheme)
    
    if (labels) {
      plot <- plot + ylab(expression(paste(Delta, Delta, Delta, "G/RT")))
    }
  }
  
  if (type == "info") {
    plot <- ggseqlogo(pwm, method = "b", font = font, col_scheme = col_scheme) +
      list(labs(x = NULL, y = NULL),
           scale_y_continuous(breaks = pretty_breaks()), alltheme)
    
    if (labels) {
      plot <- plot + ylab("Bits")
    }
  }
  
  if (type == "prob") {
    plot <- ggseqlogo(pwm, method = "c", font = font, col_scheme = col_scheme) +
      list(labs(x = NULL, y = NULL),
           scale_y_continuous(breaks = c(0.0, 0.5, 1.0), limits = c(0.0, 1.0)),
           alltheme)
    
    if (labels) {
      plot <- plot + ylab("Probability")
    }
  }
  
  suppressMessages(
    if (!axes) {
      plot <- plot + no_axes
    }
    else {
      plot <- plot + scale_x_continuous(expand = c(0.01, 0),
                                        breaks = 1:dim(pwm)[2])
    }
  )
  
  return(plot)
}

#get sequence for mesh plotting in tetrahedron
meshSeq <- function(n){
  combSeq <- expand.grid(c(1:n),c(1:n))
  combSeq <- combSeq[combSeq$Var1 > combSeq$Var2,]
  nrow(combSeq)
  seq <- c(n,n-1,n-2)
  continue <- TRUE
  while(nrow(combSeq) > 1){
    V1 <- which(combSeq$Var1 == seq[length(seq)-2])
    V2 <- which(combSeq$Var2 == seq[length(seq)-1])
    dele <- intersect(V1, V2)
    if(length(dele) > 0){
      combSeq <- combSeq[-dele,]
    }
    V1 <- which(combSeq$Var2 == seq[length(seq)-2])
    V2 <- which(combSeq$Var1 == seq[length(seq)-1])
    dele <- intersect(V1, V2)
    if(length(dele) > 0){
      combSeq <- combSeq[-dele,]
    }
    seq <- c(seq, as.numeric(combSeq[1,]))
    combSeq <- combSeq[-1,]
  }
  return(seq)
  
}

#how much is AA1>AA2 mutation in List1 explained by List2
AAChiSquare <- function(List1, List2, AA1, AA2, table = F){
  defaultW <- getOption("warn")
  options(warn = -1)
  pos1 <- which(List1 == AA1)
  pos2 <- which(List1 == AA2)
  chiSeq1 <- c()
  ChiS1 <- List1
  ChiS1[-pos1] <- '+'
  
  for(i2 in unique(List2[pos1])){
    ChiS2 <- List2
    ChiS2[ChiS2 != i2] <- '+'
    if(table){
      print(table(ChiS1, ChiS2))
    }

    chi <- chisq.test(ChiS1, ChiS2, correct = F)
    chiSeq1 <- c(chiSeq1, chi$p.value)
  }
  associatedAA1 <- unique(List2[pos1])[which(chiSeq1 == min(chiSeq1))]
  chiSeq2 <- c()
  ChiS1 <- List1
  ChiS1[-pos2] <- '+'
  for(i2 in unique(List2[pos2])){
    ChiS2 <- List2
    ChiS2[ChiS2 != i2] <- '+'
    if(table){
      print(table(ChiS1, ChiS2))
    }
    chi <- chisq.test(ChiS1, ChiS2, correct = F)
    chiSeq2 <- c(chiSeq2, chi$p.value)
  }
  associatedAA2 <- unique(List2[pos2])[which(chiSeq2 == min(chiSeq2))]
  options(warn = defaultW)
  if(associatedAA1 == associatedAA2){
    return(1 - max(c(min(chiSeq1),min(chiSeq2) )))
  }
  return(max(chiSeq))
}

#dependency between two mutations
AAmutChi <- function(List1, List2, AAfrom1, AAfrom2, AAto1, AAto2, table = F){
  defaultW <- getOption("warn")
  options(warn = -1)
  pos1 <- which(List1 == AAfrom1)
  pos2 <- which(List1 == AAto1)
  
  ChiS1 <- List1
  ChiS1[-pos1] <- '+'
  ChiS2 <- List2
  ChiS2[ChiS2 != AAfrom2] <- '+'
  if(table){
    print(table(ChiS1, ChiS2))
  }
  if(length(unique(ChiS1)) + length(unique(ChiS2)) != 4){
    chifrom <- 0.5
  }else{
    chifrom <- chisq.test(ChiS1, ChiS2, correct = F)$p.value
  }
  
  
  ChiS1 <- List1
  ChiS1[-pos2] <- '+'
  ChiS2 <- List2
  ChiS2[ChiS2 != AAto2] <- '+'
  if(table){
    print(table(ChiS1, ChiS2))
  }
  if(length(unique(ChiS1)) + length(unique(ChiS2)) != 4){
    chito <- 0.5
  }else{
    chito <- chisq.test(ChiS1, ChiS2, correct = F)$p.value
  }

  
  
  options(warn = defaultW)
  return(max(c(chifrom, chito)))
}


#form dddG table from alignment and ddG data
form.dddGList <- function(bHLH_pbAlignment, bHLH_motifs, keyPos = c(13,14,5,26,8), leaveOut = c()){
  bHLH_pbAlignment_backup <- bHLH_pbAlignment
  #start here
  bHLH_pbAlignment <- bHLH_pbAlignment_backup
  
  if(length(leaveOut) > 0){
    for(i in 1:length(leaveOut)){
      bHLH_pbAlignment <- bHLH_pbAlignment[-grep(leaveOut[i], bHLH_pbAlignment$name),]
    }
  }
  alignmentCheck <- nrow(bHLH_pbAlignment)
  if(alignmentCheck == 0){
    print('leaveOut not in Alignment')
    return()
  }
  
  bHLH_motifs_backup <- bHLH_motifs
  if(length(leaveOut > 0)){
    bHLH_motifs <- list()
    for(i in 1:length(bHLH_motifs_backup)){
      isLeave <- 0
      for(j in 1:length(leaveOut)){
        if(length(grep(leaveOut[j], bHLH_motifs_backup[[i]]$name)) > 0){
          isLeave <- 1
        }
      }
      if(isLeave == 0){
        add <- bHLH_motifs_backup[[i]]
        bHLH_motifs <- rlist::list.append(bHLH_motifs, add)
      }
    }
    if(length(bHLH_motifs) == length(bHLH_pbAlignment_backup)){
      print('nothing left out')
    }
  }
  
  #create AAtypesList
  AAtypes <- data.frame(placeHolder = c(1:nrow(bHLH_pbAlignment)))
  for(i in 1:length(keyPos)){
    AApos <- keyPos[i]
    addPos <- substr(bHLH_pbAlignment$alignment,AApos,AApos)
    AAtypes <- cbind.data.frame(AAtypes, addPos)
  }
  AAtypes <- AAtypes[,-1]
  colnames(AAtypes) <- paste0('AA', keyPos)
  
  #create dddG
  dddGList <- list()
  for(p in 1:length(keyPos)){
    AApos <- keyPos[p]
    for(posMotif in c('P1', 'P-1')){
      ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
      PosAA <- c()
      for(i in 1:length(bHLH_motifs)){
        PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
      }
      pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)
      colnames(pos_matrix) <- PosAA
      AAs <- unique(colnames(pos_matrix))
      dddGdf <- data.frame(placeHolder = c(0,0,0,0))
      for(aa1 in AAs){
        for(aa2 in AAs){
          if(aa2 != aa1){

            mean1 <- data.frame(ddG = apply(data.frame(matrix2tetrahedron(pos_matrix[,colnames(pos_matrix) == aa1])), 2, function(column) mean(column)))
            mean2 <- data.frame(ddG = apply(data.frame(matrix2tetrahedron(pos_matrix[,colnames(pos_matrix) == aa2])), 2, function(column) mean(column)))
            mean1 <- unlist(mean1)%*%ittm
            mean2 <- unlist(mean2)%*%ittm
            mean1 <- data.frame(ddG = as.numeric(mean1))
            mean2 <- data.frame(ddG = as.numeric(mean2))
        
            mean1 <- mean1 - max(mean1) + 1
            mean2 <- mean2 - max(mean2) + 1
            mean1 <- log(mean1)
            mean2 <- log(mean2)
            dddG <- data.frame(mean2 - mean1)
            
            
            colnames(dddG) <- paste0(aa1, '>', aa2)
            dddGdf <- cbind.data.frame(dddGdf, dddG)
            
          }
        }
      }
      dddGdf <- dddGdf[,-1]
      dfName <- paste0('AA', AApos, '+', posMotif)
      dddGList[[dfName]] <- dddGdf
    }
  }
  out<- list()
  out$dddGList <- dddGList
  out$AAtypes <- AAtypes
  return(out)
}

#predict energy given sequence and reference energy matrix and dddGList
predict.ddG <- function(referenceMatrix, referenceSequence, targetSequence, dddGList, posMotif = 'P-1', keyPos = c(13,14,5,8), includeChisquare = T){
  AAfrom <- c()
  AAto <- c()
  for(AApos in keyPos){
    AAfrom <- c(AAfrom,substr(referenceSequence,AApos, AApos))
    AAto <- c(AAto,substr(targetSequence,AApos, AApos))
  }
  mutDt <- data.frame(keyPos, AAfrom, AAto)
  
  deleRows <- c()
  for(i in nrow(mutDt)){
    if(mutDt$AAfrom[i] == mutDt$AAto[i]){
      deleRows <- c(deleRows, i)
    }
  }
  if(length(deleRows) == nrow(mutDt)){
    print('No difference in Key Positions')
    return(referenceMatrix)
  }else if(length(deleRows) > 0){
    mutDt <- mutDt[-deleRows,]
  }
  
  PredMatrix <- referenceMatrix
  
  for(i in 1:nrow(mutDt)){
    
    modifier <- 1
    if(i > 1 && includeChisquare){
      for(j in 1:(i-1)){
        newmod <- AAmutChi(List1 = unlist(dddGList$AAtypes[i]), List2 = unlist(dddGList$AAtypes[j]), 
                           AAfrom1 = mutDt$AAfrom[i], AAfrom2 = mutDt$AAfrom[j], AAto1 = mutDt$AAto[i], AAto2 = mutDt$AAto[j])
        if(newmod < modifier){
          modifier <- newmod
        }
      }
    }
    if(length(grep(paste0(mutDt$AAfrom[i], '>', mutDt$AAto[i]), colnames(dddGList$dddGList[[paste0('AA', mutDt$keyPos[i], '+', posMotif)]]))) != 0){
      print(modifier)
      PredMatrix <- PredMatrix + modifier * data.frame(dddGList$dddGList[[paste0('AA', mutDt$keyPos[i], '+', posMotif)]][,paste0(mutDt$AAfrom[i], '>', mutDt$AAto[i])])
    }
  }
  outMatrix <- data.frame(apply(PredMatrix, 2, function(column) column - mean(column)))
  return(outMatrix)
}

#ddG matrix of single AA 
form.ddGFeatureList <- function(bHLH_pbAlignment, bHLH_motifs, keyPos = c(13,14,5,26,8), leaveOut = c()){
  bHLH_pbAlignment_backup <- bHLH_pbAlignment
  #start here
  bHLH_pbAlignment <- bHLH_pbAlignment_backup
  
  if(length(leaveOut) > 0){
    for(i in 1:length(leaveOut)){
      bHLH_pbAlignment <- bHLH_pbAlignment[-grep(leaveOut[i], bHLH_pbAlignment$name),]
    }
  }
  alignmentCheck <- nrow(bHLH_pbAlignment)
  if(alignmentCheck == 0){
    print('leaveOut not in Alignment')
    return()
  }
  
  bHLH_motifs_backup <- bHLH_motifs
  if(length(leaveOut > 0)){
    bHLH_motifs <- list()
    for(i in 1:length(bHLH_motifs_backup)){
      isLeave <- 0
      for(j in 1:length(leaveOut)){
        if(length(grep(leaveOut[j], bHLH_motifs_backup[[i]]$name)) > 0){
          isLeave <- 1
        }
      }
      if(isLeave == 0){
        add <- bHLH_motifs_backup[[i]]
        bHLH_motifs <- rlist::list.append(bHLH_motifs, add)
      }
    }
    if(length(bHLH_motifs) == length(bHLH_pbAlignment_backup)){
      print('nothing left out')
    }
  }
  
  #create AAtypesList
  AAtypes <- data.frame(placeHolder = c(1:nrow(bHLH_pbAlignment)))
  for(i in 1:length(keyPos)){
    AApos <- keyPos[i]
    addPos <- substr(bHLH_pbAlignment$alignment,AApos,AApos)
    AAtypes <- cbind.data.frame(AAtypes, addPos)
  }
  AAtypes <- AAtypes[,-1]
  colnames(AAtypes) <- paste0('AA', keyPos)
  
  #create dddG
  dddGList <- list()
  svdList <- list()
  for(p in 1:length(keyPos)){
    AApos <- keyPos[p]
    for(posMotif in c('P1', 'P-1')){
      ProteinPos <- data.frame(name = bHLH_pbAlignment$name, pos = substr(bHLH_pbAlignment$alignment,AApos,AApos))
      PosAA <- c()
      for(i in 1:length(bHLH_motifs)){
        PosAA <- c(PosAA, ProteinPos[ProteinPos$name == bHLH_motifs[[i]]$name, 2])
      }
      pos_matrix <- gene2pos(bHLH_motifs, pos = posMotif, nrow = 4)
      colnames(pos_matrix) <- PosAA
      AAs <- unique(colnames(pos_matrix))
      dddGdf <- data.frame(placeHolder = c(0,0,0,0))
      for(aa1 in AAs){
        filteredPosMatrix <- pos_matrix[,colnames(pos_matrix) == aa1]
        mean1 <- data.frame(ddG = apply(data.frame(matrix2tetrahedron(filteredPosMatrix)), 2, function(column) mean(column)))
        mean1 <- unlist(mean1)%*%ittm
        mean1 <- data.frame(ddG = as.numeric(mean1))

        mean1 <- mean1 - max(mean1) + 1
        mean1 <- log(mean1)
        dddG <- data.frame(mean1)
        
        
        colnames(dddG) <- paste0(aa1)
        dddGdf <- cbind.data.frame(dddGdf, dddG)
        PMsvd <- svd(filteredPosMatrix)
        svdList[[paste0(AApos, '+', posMotif, '+', aa1)]] <- PMsvd
      }
      dddGdf <- dddGdf[,-1]
      dfName <- paste0('AA', AApos, '+', posMotif)
      dddGList[[dfName]] <- dddGdf
    }
  }
  out<- list()
  out$dddGList <- dddGList
  out$AAtypes <- AAtypes
  out$svd <- svdList
  return(out)
}

predEigen <- function(targetSequence, ddGFeature, keyPos = c(13,14,5,26,8), motifPos = 'P-1', print.dist = F, distType = 'Euclidean', singularity = 'one'){
  predKeyAA <- c()
  for(i in keyPos){
    predKeyAA <- c(predKeyAA, substr(targetSequence, i, i))
  }
  
  Energys <- data.frame(NULL)
  Eigens <- data.frame(NULL)
  for(i in 1:length(keyPos)){
    aa1 <- predKeyAA[i]
    pos <- keyPos[i]
    if(length(grep(aa1, colnames(ddGFeature$dddGList[[paste0('AA',pos,'+',motifPos)]]))) == 0){
      next()
    }
    Energy1 <- ddGFeature$dddGList[[paste0('AA',pos,'+',motifPos)]][,aa1]
    
    EigenM <- ddGFeature$svd[[paste0(pos, '+', motifPos, '+', aa1)]]$u
    EigenV <- ddGFeature$svd[[paste0(pos, '+', motifPos, '+', aa1)]]$d
    if(singularity == 'one'){
      Eigen1 <- EigenM[,1]
    }else if(singularity == 'all'){
      Eigen1 <- c(0,0,0,0)
      for(e in 1:ncol(EigenM)){
        Eigen1 <- Eigen1 + EigenM[,e] * EigenV[e]
      }
    }
    Energys <- rbind.data.frame(Energys, Energy1)
    Eigens <- rbind.data.frame(Eigens, Eigen1)
  }
  
  if(print.dist){
    print(dist(Eigens))
  }
  tempSum <- Energys[1,]
  for(i in 2:nrow(Energys)){
    if(distType == 'Euclidean'){
      dists <- dist(Eigens[c(1:i),])
    }else if (distType == 'Cosine'){
      dists <- c()
      cosM <- lsa::cosine(t(as.matrix(Eigens)))
      for(r in 1:nrow(cosM)){
        for(c in 1:ncol(cosM)){
          if(r > c){
            dists <- c(dists, cosM[r,c])
          }
        }
      }
      dists <- 1-dists
    }else if (distType == 'Angular'){
      dists <- c()
      cosM <- acos(lsa::cosine(t(as.matrix(Eigens))))*2/pi
      for(r in 1:nrow(cosM)){
        for(c in 1:ncol(cosM)){
          if(r > c){
            dists <- c(dists, min(cosM[r,c],1))
          }
        }
      }
    }
    
    minDist <- min(dists[(length(dists)-i+2):length(dists)])
    tempSum <- tempSum + minDist * Energys[i,]
  }
  pred <- data.frame(ddG = as.numeric(tempSum))
  PredM <- data.frame(apply(pred, 2, function(column) column - mean(column)))
  return(PredM)
}

weightedMean <- function(vector, weight){
  w <- weight/sum(weight)
  return(sum(vector*w))
}

RMSD <- function(a,b){
  return(mean((a-b)^2)^(1/2))
}


pred.SVD.ridge <- function(alignment, tetra_matrix, testSeq = '-KSLRPLLEKRRRARINQSLSQLKGLI-L------PLLGRENS--NCSKLEKADVL', keyPos = c(1:55)){
  tetra_means <- apply(tetra_matrix, 2,function(x) mean(x))
  tetra_matrix <- apply(tetra_matrix, 2,function(x) x-mean(x))
  svd <- svd(tetra_matrix)
  trueList <- c()
  predList <- c()
  #keyPos <- unique(c(X1feature,X2feature,X3feature))
  coefTable <- data.frame(NULL)
  
  #average list
  uList <- list()
  for(i in 1:3){
    rowList <- list()
    for(j in keyPos){
      u1 <- svd$u[,i]
      aa <- substr(alignment$alignment,j,j)
      dt <- cbind.data.frame(aa,u1)
      aas <- unique(dt$aa)
      means <- c()
      for(a in aas){
        subDt <- dt[dt$aa == a,]
        means <- c(means, mean(subDt$u1))
      }
      meanDt <- data.frame(aa = aas, mean = means)
      rowList[[j]] <- meanDt
    }
    uList[[i]] <- rowList
  }
  
  #synthetic U matrix
  synUList <- list()
  for(ali in 1:nrow(alignment)){
    synU <- list()
    for(i in 1:3){
      addList <- c()
      for(j in keyPos){
        aa <- substr(alignment$alignment[ali],j,j)
        uset <- uList[[i]][[j]]
        add <- uset[uset$aa == aa, 2]
        addList <- c(addList, add)
      }  
      synU[[i]] <- addList
    }
    synUList[[ali]] <- synU
  }
  
  #get linear regression coefficients
  synUpred <- list()
  coefList <- list()
  for(uindex in 1:3){
    synthesizedU <- matrix(nrow = nrow(alignment), ncol = length(keyPos))
    for(i in 1:nrow(alignment)){
      synthesizedU[i,] <- synUList[[i]][[uindex]]
    }
    synthesizedU <- data.frame(synthesizedU)
    synthesizedU$label <- svd$u[,uindex]
    cv.glm <- glmnet::cv.glmnet(as.matrix(synthesizedU[,-ncol(synthesizedU)]), synthesizedU$label, alpha = 0)
    best_lambda <- cv.glm$lambda.min
    best_model <- glmnet::glmnet(as.matrix(synthesizedU[,-ncol(synthesizedU)]), synthesizedU$label, alpha = 0, lambda = best_lambda)
    #print(best_lambda)
    summary(best_model)
    coef <- as.vector(coef(best_model))
    coef[coef == '.'] <- 0
    coefList[[uindex]] <- coef
    
    newU <- c()
    for(i in 1:nrow(alignment)){
      add <- sum(synthesizedU[i,1:length(keyPos)]*coef[1:length(keyPos)+1]) + coef[1]
      newU <- c(newU,add)
    }
    synUpred[[uindex]] <- newU
  }
  
  
  Upred <- matrix(nrow = nrow(alignment), ncol = 3)
  for(i in 1:3){
    Upred[,i] <- synUpred[[i]]
  }
  
  reSvd <- Upred %*% diag(svd$d) %*% t(svd$v)
  #print(RMSD(reSvd, tetra_matrix[-interation,]))
  
  #test
  synUtest <- list()
  for(i in 1:3){
    addList <- c()
    for(j in keyPos){
      aa <- substr(testSeq,j,j)
      uset <- uList[[i]][[j]]
      add <- uset[uset$aa == aa, 2]
      if(length(add) == 0){
        addList <- c(addList, 0)
      }else{
        addList <- c(addList, add)
      }
    }  
    synUtest[[i]] <- addList
  }
  
  predTest <- c()
  for(i in 1:3){
    predTest <- c(predTest, sum(coefList[[i]][1:length(keyPos)+1] * synUtest[[i]]) + coefList[[i]][1])
  }
  pred <- predTest %*% diag(svd$d) %*% t(svd$v)
  #print(RMSD(pred,tetra_matrix[interation,]))
  
  
  tetra_trans_matrix <- matrix(data = c(1,1,1,1,-1,-1,-1,1,-1,-1,-1,1), nrow = 4, ncol = 3, byrow = TRUE)
  pssm_trans_matrix <- MASS::ginv(tetra_trans_matrix)
  
  mean_matrix <- matrix(nrow = 1, ncol = 3, data = tetra_means, byrow = T)
  
  #apply(t((tetra_matrix+mean_matrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x)) - pos_matrix
  
  predMatrix <- matrix(nrow = 1, ncol = 3, data = pred, byrow = T)
  predMatrix1 <- apply(t((predMatrix+mean_matrix)%*%pssm_trans_matrix + 0.25), 2, function(x) x/max(x))
  
  predMatrix1[predMatrix1 < 0] <- 0.001
  
  pred <- unlist(data.frame(apply(log(predMatrix1), 2, function(column) column - mean(column))))
  pred <- as.vector(pred)
  return(pred)
}
