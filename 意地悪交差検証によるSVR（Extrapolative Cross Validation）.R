library(ClusterR) #クラスタリング
library(readr) # データ読み込み
library(dplyr) # データ操作一般
library(assertr) # データのチェック
library(rsample)　#サンプリング
library(genalg)　#遺伝的アルゴリズム
library(pls)　#PLS
library(e1071)　#SVR
library(kernlab) #マハラノビス距離計算に使用
library(iml)　#機械学習結果の可視化パッケージ
library(devtools)　#一般的各種ツール
library(parallelDist) #並列計算ツール
library(bigmemory)　#メモリ節約ツール
library(rBayesianOptimization)　#ベイズ最適化

#並列計算の設定（すでにやってある場合は実行不要）
pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
#if you want to change the number of threads for the calculation, please change the value "detectCores()"
registerDoParallel(makeCluster(detectCores()))

#BICもしくはAICに基づく最適クラスタ数の設定（クラスタリング手法はkmeansのほかにも選択可能）
ant <- as.matrix(Optimal_Clusters_KMeans(multi.regression.compounds, 15, criterion = "BIC"))
cl <- KMeans_arma(multi.regression.compounds, which(ant == ant[ant == min(ant),]), n_iter = 10, seed_mode = "random_subset",
            verbose = FALSE, CENTROIDS = NULL, seed = 1)

#各データに対し、各クラスタ中心からの距離を算出、最小値を出力する関数
minimumdistance <- function(v){
  nexa <- cl
  distances <- as.matrix(foreach(pp = 1:which(ant == ant[ant == min(ant),]), .combine = rbind,.packages = c("kernlab"))%dopar%{
    data.frame(pp, sqrt(rowSums((v - nexa[c(pp),]) ^ 2)))

  })
  distances[c(which.min(distances)), c(2)]}

#各データに対し、クラスタ中心からの最小距離＝データ密度を出力
for(i in 1:nrow(multi.regression.compounds)){
minimumdistances[i] <- minimumdistance(multi.regression.compounds[c(i), ])
}

#クラスタ中心からの最小距離が小さい順に、データベースをソート
jomoo <- minimumdistances

multi.regression.compounds <- cbind(multi.regression.compounds, jomoo)
multi.regression.compounds <- multi.regression.compounds[order(multi.regression.compounds$jomoo, decreasing = F) , ]
multi.regression.compounds <- multi.regression.compounds[, -grep(jomoo, multi.regression.compounds)]

#トレーニングとテストデータの区分（外挿側にテストデータが集まる）
train_size = 0.8

n = nrow(multi.regression.compounds)
#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds[c(1:round(n*train_size)), ]
preprocessed.y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds[-c(1:round(n*train_size)), ]
preprocessed.y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)
#---------------------------------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

#以下、SVRによるモデル構築とCVによる最適化
Xtrain <- multi.regression.x.train.s.t
ytrain <- preprocessed.y.train
datasum <- as.data.frame(cbind(ytrain,Xtrain))

xtrain = Xtrain
datasum = as.data.frame(cbind(ytrain, xtrain))

df2 <- datasum
n.validation.int <- 5
n.validation <- n.validation.int + 1
ratio.train <- 2

split.num <- n.validation + ratio.train

### SPLIT DATA INTO K FOLDS（時系列データ特論で使用したCV手法を活用、外挿側データにて検証） ###
set.seed(2016)

#cumulative = TRUEの場合、CV splitsが増えるに従って学習データも累積される
df2fold <- rolling_origin(multi.regression.compounds, initial = round(nrow(multi.regression.compounds) / split.num * ratio.train), 
                          assess = round(nrow(multi.regression.compounds) / split.num), 
                          skip = round(nrow(multi.regression.compounds) / split.num), 
                          cumulative = TRUE)
### PARAMETER LIST ###
cost <- 3
#εの候補は本来2^-10～2^0だが、今回は2^0とするとモデルが正常に構築されなかった（全て同じ値）ので0を削除した
epsilon <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1)

#グラム行列の分散を最大にするよう、γの初期値を設定
gam <- foreach(k = -20:10, .combine = rbind,.packages = c("kernlab"))%dopar%{
  rbf <- rbfdot(sigma = 2^k)
  rbf
  
  asmat <- as.matrix(xtrain)
  asmat
  
  kern <- kernelMatrix(rbf, asmat)
  sd(kern)
  data.frame((k +21), sd(kern))
}

hakata <- which.max(gam$sd.kern.)

gamma <- hakata - 21
parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {

  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  e <- parms[i, ]$epsilon
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
    deve <- df2[df2fold$splits[[j]]$in_id, ]
    test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
    mdl <- e1071::svm(ytrain~., data = deve, cost = c, gamma = 2^g, epsilon = 2^e)
    pred <- predict(mdl, test)
    data.frame(test[, c(1)], pred)
    
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
  data.frame(parms[i, ], roc)
}

epsilon <- min(result[result[, c(4)] <= (min(result[,c(4)])), c(1)])
cost <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  e <- parms[i, ]$epsilon
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
    deve <- df2[df2fold$splits[[j]]$in_id, ]
    test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
    mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
    pred <- predict(mdl, test)
    data.frame(test[, c(1)], pred)
    
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
  data.frame(parms[i, ], roc)
}

cost <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(2)])
gamma <- c(-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
### LOOP THROUGH PARAMETER VALUES ###
result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
  c <- parms[i, ]$cost
  g <- parms[i, ]$gamma
  e <- parms[i, ]$epsilon
  ### K-FOLD VALIDATION ###
  out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
    deve <- df2[df2fold$splits[[j]]$in_id, ]
    test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
    mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
    pred <- predict(mdl, test)
    data.frame(test[, c(1)], pred)
    
  }
  ### CALCULATE SVM PERFORMANCE ###
  roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
  data.frame(parms[i, ], roc)
}

gamma <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(3)])
bestperformance <- min(result[, c(4)])
bestperformance

#最適なハイパーパラメータ
gamma
cost
epsilon

#テストデータに対する予測精度検証
compounds.svr.s.t <- svm(multi.regression.x.train.s.t, preprocessed.y.train, gamma = 2^gamma, cost = 2^cost, epsilon = 2^epsilon)

ttest <- predict(compounds.svr.s.t, newdata = multi.regression.x.test.s.t)
ttrain <- predict(compounds.svr.s.t)

plot(0, 0, type = "n", xlim = c(0, 3), ylim = c(0,3),xlab = "Observed Value", ylab = "Predicted Value")

points(preprocessed.y.test, ttest, col = "orange", pch = 2)
points(preprocessed.y.train, ttrain, col = "darkgray", pch = 3)
abline(a=0, b=1)

r2ttest <- cor(preprocessed.y.test, ttest)^2
r2ttest

r2ttrain <- cor(preprocessed.y.train, ttrain)^2
r2ttrain


