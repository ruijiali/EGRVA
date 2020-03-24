library(data.table)
for (i in 1:22){
  bim=fread(paste("bcf_dup/chr",i,".mafr2.dup.bim",sep=""))
  missid=which(bim$V2==".")
  bim$V2[missid]=paste("chr",bim$V1[missid],":",bim$V4[missid],sep="")
  write.table(bim,paste("fixed_ID/chr",i,".mafr2.dup.bim",sep=""),col=F,row=F,sep="\t",quote=F)
}

rename.duplicate <- function (x, sep = "_dup", verbose = FALSE)
{
  x <- as.character(x)	
  duplix <- duplicated(x)	
  duplin <- x[duplix]	
  ix <- numeric(length = length(unique(duplin)))	
  names(ix) <- unique(duplin)	
  retval <- numeric(length = length(duplin))	
  for (i in 1:length(duplin)) {
    retval[i] <- ix[duplin[i]] <- ix[duplin[i]] + 1
  }
  
  retval <- retval + 1	
  x[duplix] <- paste(duplin, retval, sep = sep)	
  if (verbose) {
    message(sprintf("%i duplicated names", length(duplin)))
  }	
  return(list(new.x = x, duplicated.x = duplin))
}
for (i in 1:22){
  bim=fread(paste("fixed_ID/chr",i,".mafr2.dup.bim",sep=""))
  dups <- unique(bim$V2[duplicated(bim$V2)])
  ndup <- length(dups)
  if(ndup == 0)
  {
    cat("No duplicate SNPs\n")	
    q()
  } else {
    cat(paste(ndup, "duplicate SNPs\n"))
  }
  bim$V2 <- rename.duplicate(bim$V2, sep="_dup")[[1]]
  write.table(bim,paste("dup_final/chr",i,".mafr2.dup.bim",sep=""),col=F,row=F,sep="\t",quote=F)
}