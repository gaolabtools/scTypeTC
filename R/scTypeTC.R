#' predict the subtype of thyroid epithelial cell and tumor cell into thyroid follicular cell, PTC, iATC and mATC.
#' 
#' Predict subtype of single thyroid epithelial cell using multinomial-lasso model following the best vote principles.  
#' Cell with low chance of consistent prediction is labeled as undefined.   

#' @param test_data normalized gene expression data matrix from Seurat; genes in rows; cell names in columns.
#' @param cutoff cutoff on the chance of consistent predictions to report the prediction results, cell below the cutoff is predicted as undefined,default=0.5
#' @param probe_subset a set of genes as predictor in the model, automatically loaded in scTypeTC.
#' @param fit the prediction model, automatically loaded in scTypeTC.
#' @param sam.name sample name.
#' @export

scTypeTC <- function(test_data = test_data,cutoff= 0.5,sam.name="test"){
  
  ### predict: use all lambda values generated and determine the final probability with majority voting
  rownames(test_data)<- sub("-", ".",rownames(test_data)) 
  print("run prediction")
  # prepare test set
  testX = t(test_data[which(rownames(test_data) %in% rownames(probe_subset)), ])
  # fit models and predict
  pred_class = predict(fit, newx = testX, type = "class")
  # summarize probability based on all lambda values' prediction class
  iATC_prob = apply(pred_class, 1, function(x) length(which(x == "iATC"))/100)
  mATC_prob = apply(pred_class, 1, function(x) length(which(x == "mATC"))/100)
  PTC_prob = apply(pred_class, 1, function(x) length(which(x == "PTC"))/100)
  TFC_prob = apply(pred_class, 1, function(x) length(which(x == "TFC"))/100)
  pred_prob = data.frame("iATC" = iATC_prob, "mATC" = mATC_prob,
                         "PTC" = PTC_prob, "TFC" = TFC_prob)
  
  pred_class = colnames(pred_prob)[apply(pred_prob, 1, which.max)]
  pred_df = cbind(pred_prob, pred_class)
  rownames(pred_df) = rownames(testX)
  max.model <- apply(data.frame(pred_df[, 1:4]), 1,max)
  pred_confident_class <- pred_df$pred_class
  pred_confident_class[which(max.model < cutoff)] <- "undef"
  results <- data.frame(pred_df,pred_confident_class)
  results <- results[,c("TFC","PTC","iATC","mATC","pred_confident_class")]
  print("save prediction results")
  write.table(results, file=paste0(sam.name,"_final_prediction_results.txt"), row.names = TRUE, quote = FALSE, sep = "\t")
  
  print("compute shannon idex")
  smry <- table(results$pred_confident_class[which(results$pred_confident_class!="undef")])
  smry1 <- t(100*t(smry)/sum(smry))
  colnames(smry1)<-"sample"
  smry2 <- rbind(smry1[4,], smry1[3,], smry1[1,],smry1[2,])
  rownames(smry2) <- c("TFC","PTC","iATC","mATC")
  
  ##shannon index
  smry3 <- smry2/100
  sinx <- apply(smry3,2,function(x)(-sum((x+0.0001)*log((x+0.0001)))))
  print(paste0(sam.name,": shannon idex = ",sinx))
  
  
  print("make plots")
  pdf(paste0(sam.name,"_summary_barplot.pdf"), width = 5, height = 5)
  par(mfrow=c(1, 1), mar=c(5, 5, 4, 8))
  barplot(smry2,                     
          col = c("darkgreen","blue","gold2","brown"),
          legend.text = TRUE, 
          args.legend=list(
            x=ncol(smry2) + 0.7,
            y=max(colSums(smry2)),
            bty = "n"))
  dev.off()
  
  pdf(paste0(sam.name,"_summary_swirl_dotplot.pdf"), width = 5, height = 5)
  par(mfrow=c(1, 1), mar=c(2.5, 2.5, 2.5, 2.5))
    cn <- round(smry2[,1])
    colo1 <- c(rep("darkgreen",cn[1]), rep("blue",cn[2]),rep("gold2",cn[3]),rep("brown",max(0, cn[4]+100-sum(cn))))
    set.seed(001) # just to make it reproducible
    colo=sample(colo1, replace = FALSE)
    plot(1:10, rep(10,10), type="p", pch=16, ylim=c(0,11), xlim=c(0,11), cex=4, col=colo[1:10], main=colnames(smry2)[1],xaxt='n', yaxt="n",bty="n", axes = FALSE, xlab =NA,ylab=NA)
    points(1:10, rep(10,10), cex=4,lwd=2)
    points(1:10, rep(9,10), pch=16, cex=4, col=colo[11:20])
    points(1:10, rep(9,10), cex=4,lwd=2)
    points(1:10, rep(8,10), pch=16, cex=4, col=colo[21:30])
    points(1:10, rep(8,10), cex=4,lwd=2)
    points(1:10, rep(7,10), pch=16, cex=4, col=colo[31:40])
    points(1:10, rep(7,10), cex=4,lwd=2)
    points(1:10, rep(6,10), pch=16, cex=4, col=colo[41:50])
    points(1:10, rep(6,10), cex=4,lwd=2)
    points(1:10, rep(5,10), pch=16, cex=4, col=colo[51:60])
    points(1:10, rep(5,10), cex=4,lwd=2)
    points(1:10, rep(4,10), pch=16, cex=4, col=colo[61:70])
    points(1:10, rep(4,10), cex=4,lwd=2)
    points(1:10, rep(3,10), pch=16, cex=4, col=colo[71:80])
    points(1:10, rep(3,10), cex=4,lwd=2)
    points(1:10, rep(2,10), pch=16, cex=4, col=colo[81:90])
    points(1:10, rep(2,10), cex=4,lwd=2)
    points(1:10, rep(1,10), pch=16, cex=4, col=colo[91:100])
    points(1:10, rep(1,10), cex=4,lwd=2)
    
    legend("bottomright", c("ATC2","ATC1","PTC1","Normal"),
           col = c("brown","gold2","blue","darkgreen"), pch=16,
           xpd=TRUE, inset=c(0,0.9), cex=1, bty='n')
  dev.off()
  
  return(results)
  
}







