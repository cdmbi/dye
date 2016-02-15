setwd("C:/Users/Saw/Downloads/dye")

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
Substructure_fingerPrintCount <- read_file("substructure_fingerprint.csv")
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
    trainIndex <- caret::createDataPartition(data$Activity, p = .8,
                                             list = FALSE, times = 1)
    train <- data[trainIndex, ]
    test <- data[-trainIndex, ]
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




