#test for ia
miRNA = ReadmiRNA("LUAD", tumorStageList[1])
miRNA = miRNA[,order(names(miRNA))]
mRNA = ReadmRNA("LUAD", tumorStageList[1])
mRNA = mRNA[,order(names(mRNA))]

Nodes = qnorm(1- c(miRNA.LUAD.DEList[[1]]$PValue, mRNA.LUAD.DEList[[1]]$PValue))
i.names = row.names(miRNA.LUAD.DEList[[1]])
r.names = row.names(mRNA.LUAD.DEList[[1]])

Correlation = cor(miRNA[i.names[[1]],], mRNA[r.names[[1]],], method = "pearson")


Correlation = data.frame(rr = numeric(1066))
k = 0
for (i in 1:length(i.names)){
  for (j in 1:length(r.names)){
    rij = cor(miRNA[i.names[[i]],], mRNA[r.names[[j]],], method = "pearson")
    k = k + 1
    Correlation[k,] = rij
    
  }
  
}

Fr = .5*log((1+rij)/(1-rij))
