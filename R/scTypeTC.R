#' single cell sub-typing of epithelial cells in thyroid cancer with glmnet-lasso model
#'
#' Predict subtype of single thyroid epithelial cell using multinomial-lasso model following the best vote principles.
#' Cell with low chance of consistent prediction is labeled as undefined.

#' @param test_data normalized gene expression data matrix from Seurat; genes in rows; cell names in columns.
#' @param cutoff_chance cutoff on the chance of consistent predictions to report the prediction results, cell below the cutoff is predicted as undefined,default=0.5
#' @param sam.name sample name.
#' @export

#library(ComplexHeatmap)
#library(glmnet)
#library(caret)

scTypeTC <- function(test_data = test_data,cutoff_chance= 0.5,sam.name="test"){

  print("processing data")
  rownames(test_data)<- sub("-", ".",rownames(test_data))

  genes<-intersect(row.names(test_data),row.names(scdata))
  test_data<-test_data[genes,]
  scdata<-scdata[genes,]

  true.type <- true.type[order(match(true.type$cellname, colnames(scdata))), ]
  data <- data.frame(true.type$X, t(scdata))

  ##build the prediction model
  print("building the prediction model")
  ##split data
  set.seed(7)
  inTrain <- createDataPartition(y=true.type$X,p=0.3, list = FALSE)
  training <- data[inTrain,]
  training.data <- t(training[, 2:ncol(training)])
  dim(training)
  # derive subgroup information
  subgroup = training$true.type.X


  # function to run the variable selection given input data frame
  var_sel <- function(train_df, subgroup){
    # prepare data
    X = t(train_df)
    y = factor(subgroup)

    # check factor count
    summary(y)
    p = ncol(X)

    # count number of lasso non-zero coeffients
    ATC1_counts = rep(0, p)
    ATC2_counts = rep(0, p)
    PTC1_counts = rep(0, p)
    Normal_counts = rep(0, p)

    names(ATC1_counts) = colnames(X)
    names(ATC2_counts) = colnames(X)
    names(PTC1_counts) = colnames(X)
    names(Normal_counts) = colnames(X)


    alphas = seq(from=0, to = 1.0, by = 0.01)
    for (alpha in alphas) {
      message(sprintf("Run alpha = %f ...", alpha))
      # run multinomial lasso
      multinomial_lasso = glmnet(X, y, alpha = alpha, family = "multinomial", type.multinomial = "ungrouped")
      help(glmnet)
      ATC1_counts = ATC1_counts + Matrix::rowSums(multinomial_lasso$beta$ATC1 != 0)
      ATC2_counts = ATC2_counts + Matrix::rowSums(multinomial_lasso$beta$ATC2 != 0)
      PTC1_counts = PTC1_counts + Matrix::rowSums(multinomial_lasso$beta$PTC1 != 0)
      Normal_counts = Normal_counts + Matrix::rowSums(multinomial_lasso$beta$Normal != 0)
    }

    # compute selection probability
    ATC1_probs = ATC1_counts / (length(alphas) * 100)
    ATC2_probs = ATC2_counts / (length(alphas) * 100)
    PTC1_probs = PTC1_counts / (length(alphas) * 100)
    Normal_probs = Normal_counts / (length(alphas) * 100)

    probe_df = data.frame("Normal" = Normal_probs,
                          "PTC1" = PTC1_probs,
                          "ATC1" = ATC1_probs,
                          "ATC2" = ATC2_probs
    )

    return(probe_df)
  }

  training_probe = var_sel(training.data, subgroup)
  #write.table(training_probe, "1.prediction_power_of_scDEGs_30perctdata_1000runs.txt", sep="\t", quote = FALSE, row.names = TRUE)
  #training_probe <- read.delim("1.prediction_power_of_scDEGs_30perctdata_1000runs.txt", header = TRUE)
  #head(training_probe)


  # pull out any genes with prob > cut off
  prepare_input <- function(train_df, probe_subset, cutoff, test_df){

    # prepare input x y

    X = t(train_df[which(rownames(train_df) %in% rownames(probe_subset)),] )
    y = factor(subgroup)

    # prepare test set
    testX = t(test_df[which(rownames(test_df) %in% rownames(probe_subset)), ])

    return(list(X, y, testX))
  }

  prepare_input1 <- function(train_probe, cutoff){
    # subset input genes based on cutoff
    probe_subset = train_probe[which((train_probe[["Normal"]] >= cutoff) |
                                       (train_probe[["PTC1"]] >= cutoff) |
                                       (train_probe[["ATC1"]] >= cutoff) |
                                       (train_probe[["ATC2"]] >= cutoff)), ]
    #print(paste0("Cutoff=", cutoff, ", number of probes selected=", nrow(probe_subset)))

    return(probe_subset)
  }


  prepare_input3 <- function(train_probe, cutoff,low.cut){
    # subset input genes based on cutoff
    x <- apply(train_probe,1, function(x)(sort(x, decreasing = TRUE)[2]))
    probe_subset = train_probe[which(x < low.cut),]
    head(train_probe)
    n.Normal <- length(which(probe_subset$Normal>=cutoff))
    n.PTC1 <- length(which(probe_subset$PTC1>=cutoff))
    n.ATC1 <- length(which(probe_subset$ATC1>=cutoff))
    n.ATC2 <- length(which(probe_subset$ATC2>=cutoff))
    #print(paste0("low.cut=", low.cut, ", number of probes selected=", nrow(probe_subset)))
    #print(paste0("Normal predictior =", n.Normal, " genes", sep=""))
    #print(paste0("PTC1 predictior =", n.PTC1, " genes", sep=""))
    #print(paste0("ATC1 predictior =", n.ATC1, " genes", sep=""))
    #print(paste0("ATC2 predictior =", n.ATC2, " genes", sep=""))

    return(probe_subset)
  }

  cutoff<- 0.5
  low.cut <- 0.1

  temp.h=prepare_input1(training_probe, cutoff)
  temp0 =prepare_input3(temp.h,cutoff,low.cut)

  Norm.pred <- temp0[which(temp0$Normal> cutoff),]
  Norm.pred <- Norm.pred[order(Norm.pred$Normal, decreasing = TRUE),]
  PTC1.pred <- temp0[which(temp0$PTC1> cutoff),]
  PTC1.pred <- PTC1.pred[order(PTC1.pred$PTC1, decreasing = TRUE),]
  ATC1.pred <- temp0[which(temp0$ATC1> cutoff),]
  ATC1.pred <- ATC1.pred[order(ATC1.pred$ATC1, decreasing = TRUE),]
  ATC2.pred <- temp0[which(temp0$ATC2> cutoff),]
  ATC2.pred <- ATC2.pred[order(ATC2.pred$ATC2, decreasing = TRUE),]


  temp1 <- rbind(Norm.pred, PTC1.pred,ATC1.pred,ATC2.pred)
  colnames(temp1) <- c("Normal","PTC1","ATC1","ATC2")
  rownames(temp1) <- c(rownames(Norm.pred),rownames(PTC1.pred),rownames(ATC1.pred),rownames(ATC2.pred))

  ### predict: use all lambda values generated and determine the final probability with majority voting
  # print("generating the prediction results")
  predict_model2 <- function(train_df, train_probe, cutoff, test_df){
    rownames(test_df)<- sub("-", ".",rownames(test_df))
    #train_df <- training.data
    #train_probe <- temp0
    # test_df <- scdata

    # set up input matrix
    input_dfs = prepare_input(train_df, train_probe, cutoff, test_df)
    X = input_dfs[[1]]
    y = input_dfs[[2]]
    testX = input_dfs[[3]]
    dim(testX)
    # fit models
    set.seed(2021)
    fit = glmnet(X, y, family = "multinomial")
    pred_class = predict(fit, newx = testX, type = "class")

    # summarize probability based on all lambda values' prediction class
    ATC1_prob = apply(pred_class, 1, function(x) length(which(x == "ATC1"))/100)
    ATC2_prob = apply(pred_class, 1, function(x) length(which(x == "ATC2"))/100)
    PTC1_prob = apply(pred_class, 1, function(x) length(which(x == "PTC1"))/100)
    Normal_prob = apply(pred_class, 1, function(x) length(which(x == "Normal"))/100)

    pred_prob = data.frame("ATC1" = ATC1_prob, "ATC2" = ATC2_prob,
                           "PTC1" = PTC1_prob, "Normal" = Normal_prob)

    pred_class = colnames(pred_prob)[apply(pred_prob, 1, which.max)]
    #pred_class = colnames(pred_prob)[apply(pred_prob, 1, function(x)(which(x==max()))]

    pred_df = cbind(pred_prob, pred_class)
    rownames(pred_df) = rownames(testX)

    return(pred_df)
  }

  temp.h=prepare_input1(training_probe, cutoff)
  temp0 =prepare_input3(temp.h,cutoff,low.cut)


  model2.res =  predict_model2(training.data, temp0, cutoff, test_data)


  max.model <- apply(data.frame(model2.res[, 1:4]), 1,max)

  confident.50perct <- model2.res$pred_class
  confident.50perct[which(max.model < cutoff_chance)] <- "undef";

  results <- data.frame(model2.res,confident.50perct)

  print("save prediction results")
  write.table(results, file=paste0(sam.name,"_final_prediction_results.txt"), row.names = TRUE, quote = FALSE, sep = "\t")
  write.table(as.matrix(temp1),paste("prediction_power_of_DEGs_highcut_",cutoff,"_lowcut_", low.cut,".txt",sep=""),sep="\t",quote=F)

  return(results)

}







