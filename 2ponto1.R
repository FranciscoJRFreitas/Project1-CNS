#2.
#2.1)
#a)
set.seed(1111)

f <- function(x){exp(-x^2/2)/sqrt(2*pi)}

i <- integrate(f, 2, Inf)$value; i


#b)
# Integral's integrand function : f(x) = exp(-x^2/2)/sqrt(2*pi) which is similar to the normal distribution.
# P(X > 2)
p <- 1-pnorm(2,0,1); p


#c)
#Naive MC Estimator - Interval 2 to infinite.
set.seed(1111)

g <- function(y){1/y^2 * exp(-(1/y)^2/2)/sqrt(2*pi)}

m=20000
x <- runif(m,0,1/2)

I.mc1 <- (1/2-0) * mean(g(x)); I.mc1
# 0.02301548
V.mc1 <- (1/2-0)^2 * var(g(x))/m; V.mc1
# 5.579652e-08


#d)
# Naive MC Estimator - Interval 0 to 1.
set.seed(1111)

g <- function(y){((exp(-1/(2*y^2)) * (1/y^2))- (exp((-(x+1)^2)/2)))/sqrt(2*pi)}

m=20000
x <- runif(m,0,1)

I.mc2 <- (1-0) * mean(g(x)); I.mc2
# 0.02408417
V.mc2 <- var(g(x))/m; V.mc2
# 1.543204e-06

#   As x gets bigger, values tend to 0, so when we calculate the integral to the interval [0,1], we will get
# bigger results than the interval [2,Inf].


#e)
# MC Estimation
# Antithetic variable technique
set.seed(1111)

f <- function(y){((exp(-1/(2*y^2)) * (1/y^2))- (exp((-(x+1)^2)/2)))/sqrt(2*pi)}

m=20000
x=runif(m/2,0,1)

I.hat1=mean(f(x))
I.hat2=mean(f(1-x))
I.a =(I.hat1+I.hat2)/2; I.a
# 0.02287152
V.a <- 1/m*(1+cor(f(x),f(1-x)))*var(f(x)); V.a
# 4.167622e-07

# control-variate technique
set.seed(1111)

f <- function(y){((exp(-1/(2*y^2)) * (1/y^2))- (exp((-(x+1)^2)/2)))/sqrt(2*pi)}
g <- function(y){y}

m = 20000
u1 = runif(m)
cast = -(cov(u1,f(u1)))/var(u1)
cast
# -0.3850052
x = runif(m,0,1)

I.c <- mean(f(x)+cast*(g(x)-0.5)); I.c
# 0.02222363
V.c <- var(f(x))*(1-cor(f(x),g(x))^2)/m; V.c
# 1.529684e-07


#f)

100*(1-(V.c/V.mc2))
# 90.08761

100*(1-(V.a/V.mc2))
# 72.9937

res <- matrix(0,3,3)
colnames(res) <- c("Naive MC", "Ant MC", "Control MC")
rownames(res) <- c("Estimate", "Variance","% var red")
res[1,] <- round(c(I.mc2, I.a, I.c),6)
res[2,] <- round(c(V.mc2, V.a, V.c),9)
res[3,] <- c("-", round(100*(1-V.a/V.mc2),1), round(100*(1-V.c/V.mc2),1))
res

#g)
set.seed(1111);
m = 20000
h <- function(y){((exp(-1/(2*y^2)) * (1/y^2))- (exp((-(x+1)^2)/2)))/sqrt(2*pi)}
u = runif(m,0,1)
x = 2/(1-u) # 1-2/x => x = 2/(1-u)
I.hat.is <- mean(h(x)); I.hat.is
# 0.03015236
V.is <- var(h(x))/m; V.is
# 3.353631e-08

# % variance reduction to MC2
100*(1-V.is/V.mc2)
# 97.82684

# % variance reduction to MC1
100*(1-V.is/V.mc1)
# 39.89533
