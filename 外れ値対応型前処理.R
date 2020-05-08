#関連パッケージの呼び出し（必ずしもここで使うとは限りません）
library(readr) #データ読み込み
library(dplyr) #データ操作一般
library(assertr) #データのチェック
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

pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
#並列化計算で使用出来るスレッド数の設定（デフォルトは全部使用）
registerDoParallel(makeCluster(detectCores()))

#ファイルパスの読み込み（エクスプローラ出てくるので希望のcsvを選択）
path <- file.choose()
path

#CSVファイルを読み込み、compoundsという名前で保管
compounds <- read.csv(path)
View(compounds)

#不必要な列がある場合は、compounds[, -c(1,3:6)]のようにして除外（この場合は1,3～6列目を除外）
trimed.compounds <- compounds[,]

#NA値を含む列を削除（今回、1,2列目はSMILESと目的変数なのでいったん除外）
tx<-t(trimed.compounds[, -c(1,2)])
txomit<-na.omit(tx) 
trimed.compounds<-t(txomit)
View(trimed.compounds)
#compoundsの2列目（目的変数）を処理後のデータセットと結合（NA値含む行があるので、これは次の処理で削除）
trimed.compounds <- (cbind(compounds[,c(2)],trimed.compounds[, ]))

#NA値を含む行を削除
complete.compounds <- trimed.compounds[,]
is.completes <- complete.cases((complete.compounds[, ]))
is.completes
complete.compounds <- complete.compounds[is.completes,]

View(complete.compounds)

#目的変数yを設定（complete.compoundsの1列目）
y <- complete.compounds[,c(1)]
y

#外れ値を含む場合は対数変換（最小値が0となるよう）、含まない場合は正規化変換（詳細な定義は動画参照）
if (((quantile(y, 0.75) - quantile(y, 0.25))*1.5 + quantile(y, 0.75) - max(y)) < 0){                 # if ( 条件式 )
  preprocessed.y <- log((y - min(y) + 1))                   #  条件式が TRUE  のときに実行される部分
} else if ((quantile(y, 0.25) - (quantile(y, 0.75) - quantile(y, 0.25))*1.5 - min(y))*(-1) < 0){
  preprocessed.y <- (y - mean(y)) / sd(y)                   # 
}

#説明変数xを設定（complete.compoundsの1列目以外）
x <- complete.compounds[, -c(1)]

#標準偏差が0の列を削除
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]

#共相関となる列の片方を削除（相関係数の閾値はcutoffで設定、今回は0.99）
library(caret)

df2 <- as.matrix(cor(x))
hc = findCorrelation(df2, cutoff = .99, exact = FALSE) 
hc = sort(hc)
reduced_Data = x[,-c(hc)]

x <- reduced_Data
#reduced_Dataはこの後使わないので、メモリ上から削除
rm(reduced_Data)

#外れ値を含む場合は対数変換（最小値が0となるよう）、含まない場合は正規化変換（詳細な定義は動画参照）
is.greater <- apply(x, 2, function(x){((quantile(x, 0.75) - quantile(x, 0.25))*1.5 + quantile(x, 0.75) - max(x)) < 0})
is.greater

is.greater.2 <- apply(x, 2, function(x){(quantile(x, 0.25) - (quantile(x, 0.75) - quantile(x, 0.25))*1.5 - min(x))*(-1) < 0})
is.greater.2

x.g <- x[, is.greater | is.greater.2]
#x.u <- x[, -is.greater | is.greater.2]…バグフィックス（200508実施）
x.u <- x[, !(colnames(x) %in% c(colnames(x.g)))]


preprocessed.x.g <- apply(x.g, 2, function(x.g) {log((x.g - min(x.g) + 1))})
preprocessed.x.u <- apply(x.u, 2, function(x.u) {(x.u - mean(x.u)) / sd(x.u)})

preprocessed.x <- cbind(preprocessed.x.g, preprocessed.x.u)
View(preprocessed.x)

#名前の変更
multi.regression.x <- preprocessed.x[ , ]

#名前の変更及びデータフレームへの変更
multi.regression.compounds <- as.data.frame(cbind(preprocessed.y, preprocessed.x))
