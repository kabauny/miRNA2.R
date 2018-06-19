library(edgeR)

parent = "/Users/zhongningchen/Data/miRNA2"
tumorStageList = c("stage_ia", "stage_ib", "stage_iia", "stage_iib", "stage_iii", "stage_iv")

mod_filename = function(dir){
  
  df = read.table(dir, sep = "\t", header = TRUE)
  f = paste(df$id, df$filename, sep = "/")
  return (data.frame(f))
  
}
ReadmiRNA = function(cancerType, tumorStage){
  #creates dataframe given Cancer type and tumor stage. 
  #p x n
  #p name := miRNA name
  #n name := cancer_stage_index
  f = paste(parent, cancerType, sep = "/")
  
  manif = paste(f, "Manifest", paste(tumorStage, "txt", sep = "."), sep = "/")
  manifest_subset = mod_filename(manif)
  colnames(manifest_subset) = "filename"
  
  directory = paste(f, "Full", sep = "/")
  fn = paste(directory, toString(manifest_subset$filename[1]), sep = "/")
  dft = as.matrix(read.table(fn))
  dft = dft[-1,]
  
  df = data.frame(dft[,2])
  rownames(df) = dft[,1]
  colnames(df) = paste(cancerType, tumorStage, toString(1), sep = "_")
  df[,1] = as.numeric(as.character(df[,1]))
  
  for(i in 2:length(manifest_subset$filename)){
    fn = paste(directory, toString(manifest_subset$filename[i]), sep = "/")
    dft = as.matrix(read.table(fn))
    dft = dft[-1,]

    tempSampleName = paste(cancerType, tumorStage, toString(i), sep = "_")
    df[tempSampleName] = as.numeric(as.character(dft[,2]))
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
  
  qlf = glmQLFTest(fit, coef = 2)
  return(qlf$table[qlf$table$PValue < 0.05,])
  
  #tr = glmTreat(fit, coef = 2, lfc = 1)
  #return(tr$table[tr$table$PValue < 0.05,])
}
stage.compare2 = function(comp, NORM, stage, ts){
  comp = comp[stage$state == ts] 
  NORM = NORM[stage$state == ts]
  filter = Zero_filter(comp, NORM)
  comp.filter = comp[filter,]
  NORM.filter = NORM[filter,]
  
  return(DE(comp.filter, NORM.filter))
  
}
#stage.compare3 = function(comp, NORM, stage){

  filter = Zero_filter(comp, NORM)
  comp.filter = comp[filter,]
  NORM.filter = NORM[filter,]
  
  return(stage.comparison(comp.filter, NORM.filter))
  
}

miRNA.LUAD.comp = ReadmiRNA("LUAD", "complement")
miRNA.LUAD.normal = ReadmiRNA("LUAD", "normal")
stage.LUAD = compStage("LUAD")

filter = Zero_filter(miRNA.LUAD.comp, miRNA.LUAD.normal)
miRNA.LUAD.comp.filtered = miRNA.LUAD.comp[filter,]
miRNA.LUAD.normal.filtered = miRNA.LUAD.normal[filter,]
dfexp.total = DE(df_temp = miRNA.LUAD.comp.filtered, normal = miRNA.LUAD.normal.filtered)

ia.DE = stage.compare2(miRNA.LUAD.comp, miRNA.LUAD.normal, stage.LUAD, tumorStageList[1])
ib.DE = stage.compare2(miRNA.LUAD.comp, miRNA.LUAD.normal, stage.LUAD, tumorStageList[2])
iia.DE = stage.compare2(miRNA.LUAD.comp, miRNA.LUAD.normal, stage.LUAD, tumorStageList[3])
iib.DE = stage.compare2(miRNA.LUAD.comp, miRNA.LUAD.normal, stage.LUAD, tumorStageList[4])
iii.DE = stage.compare2(miRNA.LUAD.comp, miRNA.LUAD.normal, stage.LUAD, tumorStageList[5])

sum(stage.LUAD$state == tumorStageList[1])
sum(stage.LUAD$state == tumorStageList[2])
sum(stage.LUAD$state == tumorStageList[3])
sum(stage.LUAD$state == tumorStageList[4])
sum(stage.LUAD$state == tumorStageList[5])
