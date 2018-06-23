library(edgeR)

parent = "/Users/zhongningchen/Data/mRNA"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")

mod_filename = function(dir){
  
  df = read.table(dir, sep = "\t", header = TRUE)
  f = paste(df$id, df$filename, sep = "/")
  return (data.frame(f))
  
}
ReadmRNA = function(cancerType, tumorStage){
  #creates dataframe given Cancer type and tumor stage. 
  #p x n
  #p name := miRNA name
  #n name := cancer_stage_index
  parent = "/Users/zhongningchen/Data/mRNA"
  f = paste(parent, cancerType, sep = "/")
  
  manif = paste(f, "Manifest", paste(tumorStage, "txt", sep = "."), sep = "/")
  manifest_subset = mod_filename(manif)
  colnames(manifest_subset) = "filename"
  
  directory = paste(f, "Full", sep = "/")
  fn = paste(directory, toString(manifest_subset$filename[1]), sep = "/")
  dft = read.table(gzfile(fn))
  
  dft = head(dft, -5)
  
  df = data.frame(dft$V2)
  rownames(df) = dft$V1
  
  f2 = paste(f, "Manifest", paste("one2one", "txt", sep = "."), sep = "/")
  df2 = as.matrix(read.table(f2))[-1,] #df2[,3] is miRNA, df2[,4] is pt ID
  w = which(df2[,2] == toString(manifest_subset$filename[1]))
  cname = df2[w, 4]
  
  colnames(df) = cname
  
  
  for(i in 2:length(manifest_subset$filename)){
    fn = paste(directory, toString(manifest_subset$filename[i]), sep = "/")
    dft = read.table(gzfile(fn))
    dft = head(dft, -5)
    
    
    if (length(which(df2[,2] == toString(manifest_subset$filename[i]))) > 0){
      w = which(df2[,2] == toString(manifest_subset$filename[i]))
      tempSampleName = df2[w, 4]
      df[tempSampleName] = as.numeric(as.character(dft[,2]))
      
    }

  }
  
  return(df)
}
compStage = function(cancerType){
  f = paste(parent, cancerType, sep = "/")
  manif = paste(f, "Manifest", paste("complement", "txt", sep = "."), sep = "/")
  df = read.table(manif, sep = "\t", header = TRUE)
  return(data.frame(df["state"]))
  
}
Zero_filter = function(df1, df2){
  comp1 = apply(df1, 1, function(x) all(x > 10))
  comp2 = apply(df2, 1, function(x) all(x > 10))
  comp3 = (comp1 | comp2)
  
  return(comp3)
  
}
DE = function(df_temp, normal){
  
  tmp1 = normal
  tmp2 = df_temp
  
  x = cbind(tmp1, tmp2)
  group = c(rep(1, dim(tmp1)[2]), rep(2, dim(tmp2)[2]))
  y = DGEList(counts = x, group = group, genes = rownames(x))
  o = order(rowSums(y$counts), decreasing = TRUE)
  y = y[o,]
  y = calcNormFactors(y)
  patient = factor(c(1:dim(tmp1)[2], 1:dim(tmp2)[2]))
  tissue = factor(c(rep("Normal", dim(tmp1)[2]), rep("LUAD", dim(tmp2)[2])))
  design = model.matrix(~patient+tissue)
  y = estimateDisp(y, design, robust= TRUE)
  
  #classic: 
  
  #fit = glmFit(y, design)
  #lrt = glmLRT(fit)
  #return(topTags(lrt, n = Inf, p.value = 1e-5))
  
  fit = glmQLFit(y, design)
  
  #qlf = glmQLFTest(fit, coef = 2)
  #return(qlf$table[qlf$table$PValue < 0.05,])
  
  tr = glmTreat(fit, coef = 2, lfc = 1)
  return(tr$table[tr$table$PValue < 0.05,])
}
stage.compare2 = function(comp, NORM, stage, ts){
  comp = comp[stage$state == ts] 
  NORM = NORM[stage$state == ts]
  filter = Zero_filter(comp, NORM)
  comp.filter = comp[filter,]
  NORM.filter = NORM[filter,]
  
  return(DE(comp.filter, NORM.filter))
  
}

mRNA.LUAD.comp = ReadmRNA("LUAD", "complement")
mRNA.LUAD.normal = ReadmRNA("LUAD", "normal")[,-c(57:59)]
mRNA.stage.LUAD = compStage("LUAD")
filter = Zero_filter(mRNA.LUAD.comp, mRNA.LUAD.normal)
mRNA.LUAD.comp.filtered = mRNA.LUAD.comp[filter,]
mRNA.LUAD.normal.filtered = mRNA.LUAD.normal[filter,]
mRNA.DE = DE(df_temp = mRNA.LUAD.comp.filtered, normal = mRNA.LUAD.normal.filtered)


#differential expression of mRNA 
mRNA.ia.DE = stage.compare2(mRNA.LUAD.comp, mRNA.LUAD.normal, mRNA.stage.LUAD, tumorStageList[1])
mRNA.ib.DE = stage.compare2(mRNA.LUAD.comp, mRNA.LUAD.normal, mRNA.stage.LUAD, tumorStageList[2])
mRNA.iia.DE = stage.compare2(mRNA.LUAD.comp, mRNA.LUAD.normal, mRNA.stage.LUAD, tumorStageList[3])
mRNA.iib.DE = stage.compare2(mRNA.LUAD.comp, mRNA.LUAD.normal, mRNA.stage.LUAD, tumorStageList[4])
mRNA.iii.DE = stage.compare2(mRNA.LUAD.comp, mRNA.LUAD.normal, mRNA.stage.LUAD, tumorStageList[5])

#differential expression ordered by P value
mRNA.ia.DE = mRNA.ia.DE[order(mRNA.ia.DE$PValue),]
mRNA.ib.DE = mRNA.ib.DE[order(mRNA.ib.DE$PValue),]
mRNA.iia.DE = mRNA.iia.DE[order(mRNA.iia.DE$PValue),]
mRNA.iib.DE = mRNA.iib.DE[order(mRNA.iib.DE$PValue),]
mRNA.iii.DE = mRNA.iii.DE[order(mRNA.iii.DE$PValue),]

names_only.mRNA1 = rownames(mRNA.ia.DE)
names_only.mRNA2 = rownames(mRNA.ib.DE)
names_only.mRNA3 = rownames(mRNA.iia.DE)
names_only.mRNA4 = rownames(mRNA.iib.DE)
names_only.mRNA5 = rownames(mRNA.iii.DE)




# names_only.mRNA.ia.DE = gsub("\\..*", "",rownames(mRNA.ia.DE))
# names_only.mRNA.ib.DE = gsub("\\..*", "",rownames(mRNA.ib.DE))
# names_only.mRNA.iia.DE = gsub("\\..*", "",rownames(mRNA.iia.DE))
# names_only.mRNA.iib.DE = gsub("\\..*", "",rownames(mRNA.iib.DE))
# names_only.mRNA.iii.DE = gsub("\\..*", "",rownames(mRNA.iii.DE))

write.csv(names_only.mRNA.ia.DE, "names_only_mRNA_stage_ia_DE.csv", row.names = F)
write.csv(names_only.mRNA.iii.DE, "names_only_mRNA_stage_iii_DE.csv", row.names = F)



