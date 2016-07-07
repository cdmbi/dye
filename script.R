#setwd("C:/Users/Saw/Downloads/dye")
#### function for reading file
read_file <- function(x){
  library(caret)
  library(data.table)
  data <- fread(x)
  data <- as.data.frame(data)
  descriptors <- data[, 2:ncol(data)]
  set.seed(1)
  filtered_descriptors <- descriptors[, -nearZeroVar(descriptors)]
  Activity <- data$Activity
  filtered_data <- cbind(Activity, filtered_descriptors)
  return(filtered_data)
}



AtomPairs2D_fingerPrintCount <- read_file("AtomPairs2DFingerprintCount.csv")
AtomPairs2D_fingerPrinter <- read_file("AtomPairs2DFingerprinter.csv")
Substructure_fingerPrintCount <- read_file("substructure_fingerprint_count.csv")
Substructure_fingerPrinter <- read_file("substructure_fingerprint.csv")
Extended_finterPrinter <- read_file("extended_finger_printer.csv")
FingerPrinter <- read_file("finger_printer.csv")
Estate_FingerPrinter <- read_file("estate_fingerprint.csv")
GraphOnly_FingerPrinter <- read_file("graph_only_fingerprint.csv")
KlekotaRoth_FingerprintCount <- read_file("KlekotaRothFingerprintCount.csv")
KlekotaRoth_FingerPrinter <- read_file("KlekotaRothFingerprinter.csv")
MACCS_FingerPrinter <- read_file("MACCS_finger_printer.csv")
Pubchem_FingerPrinter <- read_file("Pubchem_finger_printer.csv")


input <- list(AtomPairs2D_fingerPrintCount=AtomPairs2D_fingerPrintCount,
              AtomPairs2D_fingerPrinter = AtomPairs2D_fingerPrinter,
              Substructure_fingerPrintCount = Substructure_fingerPrintCount,
              Substructure_fingerPrinter = Substructure_fingerPrinter,
              Extended_finterPrinter = Extended_finterPrinter,
              FingerPrinter = FingerPrinter,
              Estate_FingerPrinter = Estate_FingerPrinter,
              GraphOnly_FingerPrinter = GraphOnly_FingerPrinter,
              KlekotaRoth_FingerprintCount = KlekotaRoth_FingerprintCount,
              KlekotaRoth_FingerPrinter = KlekotaRoth_FingerPrinter,
              MACCS_FingerPrinter = MACCS_FingerPrinter,
              Pubchem_FingerPrinter = Pubchem_FingerPrinter)


RF_training_classification <- function(x) {
  library(parallel)
  library(doSNOW)
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  
  ok <- list(100)
  ok <- foreach(i = 1:100 ) %dopar% {
    
    data <- x
    #trainIndex <- caret::createDataPartition(data$Activity, p = .8,
    #                                         list = FALSE, times = 1)
    #train <- data[trainIndex, ]
    #test <- data[-trainIndex, ]
    train <- data
    model_train <- ranger::ranger(Activity~., data = train, write.forest = TRUE, save.memory = TRUE)
    rm(ctrl)
    rm(rf)
    rm(data)
    rm(trainIndex)
    actual <- train$Activity
    prediction <- predict(model_train, train)
    prediction <- prediction$predictions
    rm(train)
    rm(test)
    rm(model_train)
    results <- caret::confusionMatrix(prediction, actual)
    results <- results$table
    results <- as.numeric(results)
    rm(prediction)
    rm(actual)
    ok[[i]] <- results
    
    
    
  }
  return(ok)
  stopCluster(cl)
}


mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 2),
    round(sd(x, na.rm = TRUE), digits = 2))
}

RF_train_classification <- function(x) {
  ok <- RF_training_classification(x)
  results <- data.frame(ok)
  rm(ok)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  rm(m)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}


#training <- randomForest_train(data)

training_results_classification <- lapply(input, function(x) {
  results <- RF_train_classification(x)
  return(results)
})

RF_10_CV <- function(x){
  library(parallel)
  library(doSNOW)
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  
  ok <- list(100)
  ok <- foreach(i = 1:100, .packages = 'randomForest') %dopar% { 
    data <- x
    #trainIndex <- caret::createDataPartition(data$Activity, p = .8,
     #                                        list = FALSE, times = 1)
   # train_control <- caret::trainControl(method="LOOCV")
    # train the model
    model <- randomForest::randomForest(Activity~., data = data)
    #model <- caret::train(Activity~., data=data, trControl=train_control, method="rf")
    prediction <- predict(model, data)
    actual <- data$Activity
    
    results <- caret::confusionMatrix(prediction, actual)
    results <- results$table
    results <- as.numeric(results)
    rm(data)
    rm(model)
    rm(prediction)
    rm(myData)
    rm(actual)
    ok[[i]] <- results
  }
  return(ok)
  stopCluster(cl)
}


mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 2),
    round(sd(x, na.rm = TRUE), digits = 2))
}

RF_10_cross_validation <- function(x) {
  ok <- RF_10_CV(x)
  results <- data.frame(ok)
  rm(ok)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  return(results_all)
}


cross_validation_results_classification <- lapply(input, function(x) {
  results <- RF_10_cross_validation(x)
  return(results)
})


RF_training_classification_external <- function(x) {
  library(parallel)
  library(doSNOW)
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  
  ok <- list(100)
  ok <- foreach(i = 1:100 ) %dopar% {
    
    data <- x
    trainIndex <- caret::createDataPartition(data$Activity, p = .8,
                                             list = FALSE, times = 1)
    train <- data[trainIndex, ]
    test <- data[-trainIndex, ]
    model_train <- ranger::ranger(Activity~., data = train, write.forest = TRUE, save.memory = TRUE)
    rm(ctrl)
    rm(rf)
    actual <- test$Activity
    prediction <- predict(model_train, test)
    prediction <- prediction$predictions
    rm(train)
    rm(model_train)
    results <- caret::confusionMatrix(prediction, actual)
    results <- results$table
    results <- as.numeric(results)
    rm(prediction)
    rm(actual)
    ok[[i]] <- results
    
    
    
  }
  return(ok)
  stopCluster(cl)
}


mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 2),
    round(sd(x, na.rm = TRUE), digits = 2))
}

RF_train_classification <- function(x) {
  ok <- RF_training_classification_external(x)
  results <- data.frame(ok)
  rm(ok)
  data <- data.frame(results)
  rm(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  rm(ACC)
  rm(SENS)
  rm(SPEC)
  rm(MCC)
  rm(data)
  rm(m)
  results_all <- (data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
  rownames(results_all) <- c("ACC_Mean", "ACC_SD", "Sens_Mean", "Sens_SD", "Spec_Mean", "Spec_SD",
                             "MCC_Mean", "MCC_SD")
  rm(results_ACC)
  rm(results_SENS)
  rm(results_SPEC)
  rm(results_MCC)
  return(results_all)
}


#testing <- randomForest_train(data)

testing_results_classification <- lapply(input, function(x) {
  results <- RF_train_classification(x)
  return(results)
})


randomForest_feature_importance <- function(x) {
  library(doSNOW)
  library(foreach)
  library(parallel)
  cl <- makeCluster(8)
  registerDoSNOW(cl)
  
  results <- list(100)
  results <- foreach (i = 1:100) %dopar% {
    
    data <- x
    trainIndex <- caret::createDataPartition(data$Activity, p = .9,
                                             list = FALSE, times = 1)
    train <- data[trainIndex, ]
    test <- data[-trainIndex, ]
    model <- randomForest::randomForest(Activity~., data = train, importance = TRUE)
    #sel <- prospectr::naes(x, k = 90, pc = 5, iter.max = 100)
    #myData <- x[sel$model, ]
    #Test <- x[sel$test, ]
    #ctrl <- caret::trainControl(method = "repeatedcv", number = 10, repeats = 1)
    #tune <- caret::train(pIC50 ~., data = myData,  method = "rf", tuneLength = 10,
    #                     trControl = ctrl)
    #model <- randomForest::randomForest(pIC50~., data = myData, importance = TRUE,
    #                                    mtry = tune$bestTune[[1]])
    importance <- model$importance
    rm(data)
    rm(trainIndex)
    rm(trian)
    rm(test)
    rm(model)
    results[[i]] <- importance
  }
  return(results)
  stopCluster(cl)
}

plot_importance <- function(x) {
  library(ggplot2)
  library(cowplot)
  randomForest_feature_importance_result <- randomForest_feature_importance(x)
  df <- as.data.frame(randomForest_feature_importance_result)
  df <- t(apply(df, 1, round, digits = 4))
  index <- seq(2, 200, by = 2)
  df <- df[, index]
  mean <- data.frame(apply(df, 1, mean))
  sd <- data.frame(apply(df, 1, sd))
  geni_index <- cbind(mean, sd)
  colnames(geni_index) <- c("mean", "sd")
  geni_index <- geni_index[order(mean, decreasing = TRUE),]
  geni_index <- head(geni_index, 30)
  set.seed(200)
  a <- cbind(descriptors = rownames(geni_index), geni_index)
  a$descriptors <- factor(a$descriptors, levels = a[order(a$mean), "descriptors"])
  z <- ggplot2::ggplot(a, aes(x = mean, y = descriptors)) +
    geom_point(size = 4, colour = "black", fill = "red", pch = 21) +
    geom_errorbarh(aes(xmin = mean-sd, xmax = mean+sd), colour = "black") +
    ggtitle("") + xlab("Gini index") + ylab ("") +
    theme(
      plot.margin = grid::unit(c(1, 1, 1, 1), "cm"),
      panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
      axis.title.x = element_text(size = 20, face = "bold", colour = "black"))
  z + cowplot::background_grid(major = "xy", minor = "none")
      
    
}

save_pdf <- function(x) {
  
}

pdf("sampleGraph", width = 5, height = 10)
dev.off()



a <- plot_importance(AtomPairs2D_fingerPrinter)
b <- plot_importance(AtomPairs2D_fingerPrintCount)
c <- plot_importance(Substructure_fingerPrinter)
d <- plot_importance(Substructure_fingerPrintCount)
e <- plot_importance(Extended_finterPrinter)
f <- plot_importance(FingerPrinter) #### re start
g <- plot_importance(Estate_FingerPrinter)
h <- plot_importance(GraphOnly_FingerPrinter)
i <- plot_importance(KlekotaRoth_FingerprintCount)
j <- plot_importance(KlekotaRoth_FingerPrinter)
k <- plot_importance(MACCS_FingerPrinter)
l <- plot_importance(Pubchem_FingerPrinter)



Extended_finterPrinter <- read_file("extended_finger_printer.csv")
FingerPrinter <- read_file("finger_printer.csv")
Estate_FingerPrinter <- read_file("estate_fingerprint.csv")
GraphOnly_FingerPrinter <- read_file("graph_only_fingerprint.csv")
KlekotaRoth_FingerprintCount <- read_file("KlekotaRothFingerprintCount.csv")
KlekotaRoth_FingerPrinter <- read_file("KlekotaRothFingerprinter.csv")
MACCS_FingerPrinter <- read_file("MACCS_finger_printer.csv")
Pubchem_FingerPrinter <- read_file("Pubchem_finger_printer.csv")



















geni_index <- geni_index[order(mean, decreasing = TRUE),]
geni_index <- head(geni_index, 20)
set.seed(200)
a <- cbind(descriptors = rownames(geni_index), geni_index)
a$descriptors <- factor(a$descriptors, levels = a[order(a$mean), "descriptors"])
### plot
z <- ggplot(a, aes(x = mean, y = descriptors)) +
  geom_point(size = 4, colour = "black", fill = "red", pch = 21) +
  geom_errorbarh(aes(xmin = mean-sd, xmax = mean+sd), colour = "black") +
  ggtitle("") + xlab("Gini index") + ylab("") +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    panel.border = element_rect(linetype = "solid", colour = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 20, face = "bold", colour = "black"))
z + background_grid(major = "xy", minor = "none")

