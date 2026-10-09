library(ggplot2)
library(cowplot)
set.seed(1)
n <- 10000
a <- 10 # shape
s <- 1/a # scale
# s <- 2000 (for a = 0.1)
x <- rgamma(n,shape = a,scale = s)
print(mean(x))
print(median(x))
print(qgamma(0.5,shape = a,scale = s))
print(sd(x))
hist(log2(x),n = 64)
stop()
par(mfcol = c(1,2))
print(sd(x))
s <- 1.15
x <- rnorm(n,0,s)
y <- exp(x)
hist(y,n = 64)
print(sd(y))
y <- seq(-2,2,length.out = 1000)
x <- 10^y
x1 <- rgamma(n,shape = 0.1,scale = 2000)
x2 <- rgamma(n,shape = 10,scale = 0.1)
par(mfcol = c(2,1))
hist(log2(x1),n = 64,col = "red",xlim = c(-20,20))
hist(log2(x2),n = 64,col = "blue",xlim = c(-20,20))

ggplot(diamonds,aes(price, fill = cut)) +
  geom_histogram(binwidth = 500,position = "stack")
ggplot(diamonds,aes(price, fill = cut)) +
  geom_histogram(binwidth = 500,position = "identity")
