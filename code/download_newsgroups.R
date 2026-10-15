# Copied newsgroups prep code from https://stephenslab.github.io/fastTopics-experiments/prepare_newsgroups.html
# Downloaded matlab version gz & vocabulary.txt
library(Matrix)
library(readr)
library(tools)
data.dir <- here("docs","data","20news-bydate")

topic.names   <- read.table(file.path(data.dir,"train.map"),sep = " ",
                            stringsAsFactors = FALSE)[[1]]
words         <- read.table(file.path(data.dir,"vocabulary.txt"),
                            stringsAsFactors = FALSE)[[1]]
topics        <- c(read.table(file.path(data.dir,"train.label"))[[1]],
                   read.table(file.path(data.dir,"test.label"))[[1]])
train         <- read.table(file.path(data.dir,"train.data"),sep = " ")
test          <- read.table(file.path(data.dir,"test.data"),sep = " ")
names(train)  <- c("document","word","count")
names(test)   <- c("document","word","count")
test$document <- test$document + max(train$document)
dat           <- rbind(train,test)

counts           <- sparseMatrix(i = dat$document,j = dat$word,x = dat$count)
colnames(counts) <- words

cols   <- which(colSums(counts > 0) >= 2)
counts <- counts[,cols]

topics         <- factor(topics)
levels(topics) <- topic.names

n <- nrow(counts)
m <- ncol(counts)
cat(sprintf("Number of documents: %d\n",n))
cat(sprintf("Number of terms in vocabulary: %d\n",m))
cat(sprintf("Rate of nonzero counts: %0.2f%%\n",100*nnzero(counts)/(n*m)))

cat("The word counts are mostly small, with a small number of large counts:\n")
print(quantile(summary(counts)$x,c(0,0.5,0.9,0.99,0.999,1)))

save(list = c("counts","topics"),file = "newsgroups.RData")
resaveRdaFiles("newsgroups.RData")
