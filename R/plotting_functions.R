



plot_gene_bySubtype <- function(data,dataName,subtypeMethod,gene){
  par(mai=c(1,1,1,1))

  colnames(data) <- c("subtypeMethod","geneName")
  
  boxplot(geneName ~ subtypeMethod,data=data,main=paste(gene,"expression across",dataName,"samples","-",subtypeMethod,", N:",dim(data)[1])
          ,ylab = paste(gene,"mRNA level\n log2(TPM+0.001)"), col=c("#de2d26",rev(c('#eff3ff','#bdd7e7','#6baed6','#3182bd')))
          ,outline=F, ylim=c(min(data$geneName,na.rm = T),max(data$geneName,na.rm = T)))
  
  stripchart(geneName ~ subtypeMethod, vertical = TRUE, data = data, 
             method = "jitter", add = TRUE, pch = 20, col=rgb(0, 0, 0, 0.4),jitter=0.2)
  
  tests <- names(table(data$subtypeMethod))
  
  pvals <- c()
  for (i in 1:(length(tests)-1) ){
    for (j in (i+1):length(tests)) {
      pvals <- c(pvals,kruskal.test(list(data[data$subtypeMethod==tests[i],"geneName"],data[data$subtypeMethod==tests[j],"geneName"]))$p.value)
    }
  }
  
  Colpvals <- c()
  for (i in 1:(length(tests)-1) ){
    for (j in (i+1):length(tests)) {
      Colpvals <- c(Colpvals,paste(tests[i],tests[j],sep = "_"))
    }
  }
  names(pvals) <- Colpvals
  
  numbers <- unlist(lapply(tests, function(test){
    sum(data$subtypeMethod==test)
  }))
  
  names(numbers) <- tests
  
  return(list(pvals,numbers))
  
}

plot_gene_bySubtype_specific <- function(data,dataName,subtypeMethod,subtypeOfInterest,gene){

  colnames(data) <- c("subtypeMethod","geneName")
  data[which(data$subtypeMethod==subtypeOfInterest),"groups"] <- subtypeOfInterest
  data[which(data$subtypeMethod!=subtypeOfInterest),"groups"] <- paste("non-",subtypeOfInterest,sep = "")
  print(table(data[,"groups"]))
  
  
  A <- kruskal.test(list(data[which(data$subtypeMethod==subtypeOfInterest),"geneName"],data[which(data$subtypeMethod!=subtypeOfInterest),"geneName"]))
  # wilcox.test(x=data[which(data$subtypeMethod==subtypeOfInterest),"geneName"],y=data[which(data$subtypeMethod!=subtypeOfInterest),"geneName"],exact = T,correct = T)
  boxplot(geneName ~ groups,data=data[,c("geneName","groups")],main=paste(names(gene),"expression across",dataName,"samples","-",subtypeMethod
                                                                          ,"\nN:",dim(data)[1],paste(", pval:",sprintf("%.2E",A$p.value))
                                                                          ,",",names(table(data[,"groups"]))[1],":",table(data[,"groups"])[1]
                                                                          ,",",names(table(data[,"groups"]))[2],":",table(data[,"groups"])[2])
          ,outline=F
          , ylim=c(min(data$geneName,na.rm = T),max(data$geneName,na.rm = T))
          ,ylab = paste(gene,"mRNA level\nlog2(TPM+0.001)"),col=c("#de2d26",'#2171b5'))
  stripchart(geneName ~ groups, vertical = TRUE, data = data[,c("geneName","groups")], 
             method = "jitter", add = TRUE, pch = 20, col=rgb(0, 0, 0, 0.4),jitter=0.2)
  #legend("topright",legend = paste("pval:",sprintf("%.2E",A$p.value)),bty = "n")
  
}
