##1.1

## (b)
pdf.pareto = function(x, alpha, L,H){
    (alpha * L^(alpha) * x^(-alpha-1))/(1 - (L/H)^(alpha))
}

sim.IT = function(m,alpha, L, H){
    if(L>H) stop("Lower bound greater than Higher bound")
    else if(alpha <= 0) stop("Alpha must be a positive number")
    else if(L <= 0) stop("Lower Bound must be a positive number")
    else if(m <= 0) stop("m must be a positive number")
    x=0
    t=0
    while(t<m){
        t=t+1
        u=runif(1,0,1)
        up = -(u*(1-((L/H)^alpha))-1)
        down = L^alpha
        x[t]= (up/down)^(1/(-alpha))
        
    }
    x    
}

## (d)
library(ggplot2)
set.seed(2447)
alpha = 0.25
L = 2
H = 3
m = 12000
set.IT = sim.IT(m,alpha,L,H)
pareto <- function(x){pdf.pareto(x,alpha,L,H)}
p1 <- ggplot(data.frame(sample = set.IT[1:15]), aes(x = sample)) +
     geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.2,
                   breaks=seq(L,H,0.1)) +
                   geom_function(fun=pareto)
plot(p1)

## (e)
sim.AR = function(m,alpha, L,H){
    if(L>H) stop("Lower bound greater than Higher bound")
    else if(alpha <= 0) stop("Alpha must be a positive number")
    else if(L <= 0) stop("Lower Bound must be a positive number")
    else if(m <= 0) stop("m must be a positive number")
    x = vector()
    x_rej = vector()
    y = vector()
    y_rej = vector()
    f =function (x){ (alpha * L^(alpha) * x^(-alpha-1))/(1 - (L/H)^(alpha)) }
    g = function(x){alpha*exp((-alpha)*x)}
    h = function(x){f(x)/g(x)}
    M = optimize(h, c(L,H),maximum = T)$objective
    g.CDF=function(x){1-exp(-alpha*x)}
    g.ITM = function(x){log(1-x)/(-alpha)}
    for(i in 1:m){
        u <- 1
        a <- 0
        while(u> a){
            x.c   <- g.ITM(runif(1,g.CDF(2),g.CDF(3))) # O g.CDF aqui serve
            # para não estar a gerar valores desnecessários 
            # através da distribuição exponencial
            a <- f(x.c)/(M*g(x.c)) 
            u <- runif(1,0,1)
            if(u > a) {
                x_rej = c(x_rej, x.c )
                y_rej = c(y_rej, u*M*g(x.c))
            }
        }
        x <- c(x,x.c)
        y <- c(y,u*M*g(x.c))
    }
    list(x = x,x_rej = x_rej,y = y, y_rej =y_rej,
         failRate= length(x_rej)/(length(x_rej) + length(x)))
}

## (g)

alpha = 0.25
L = 2
H = 3
m=12000
f =function (x){ (alpha * L^(alpha) * x^(-alpha-1))/(1 - (L/H)^(alpha)) }
g = function(x){alpha*exp((-alpha)*x)}
h = function(x){f(x)/g(x)}
M = optimize(h, c(L,H),maximum = T)$objective
Mg= function(x){M*g(x)}
set.AR = sim.AR(m,alpha,L,H)
pareto <- function(x){pdf.pareto(x,alpha,L,H)}

library(ggplot2)
library(gridExtra)
p2 <- ggplot(data.frame(sample = set.AR$x[1:15]), aes(x = sample)) +
     geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.2,
                   breaks=seq(L,H,0.1)) +
    geom_function(fun=pareto)
p3 <- ggplot( ) +
    geom_point(data.frame(list(x = set.AR$x,y =set.AR$y)),
               mapping = aes(x = x, y=y),
               size = 0.05,
               colour = "black") +
    geom_point(data.frame(list(x_rej = set.AR$x_rej,y_rej =set.AR$y_rej)),
               mapping = aes(x = x_rej, y=y_rej),
               size = 0.05,
               colour = "red") +
    xlim(L,H) + geom_function(fun=f, aes(colour = "f(x)")) +
    geom_function(fun=Mg, aes(colour = "Mg(x)")) 

grid.arrange(p2, p3, ncol = 2)

## (h)

compareTimes = function(m,nruns,alpha,L,H){
    library(Pareto)
    v.IT = vector()
    v.AR = vector()
    v.Pareto = vector()
    i=0
    while(i < nruns){
#        print(i)
        ptm <- proc.time()
        z = sim.IT(m,alpha,L,H)
        v.IT = c(v.IT,proc.time() -ptm)
        ## ptm <- proc.time()
        ## z = sim.AR(m,alpha,L,H)
        ## v.AR = c(v.AR,proc.time() -ptm)
        ptm <- proc.time()
        z=rPareto(m,L,alpha,truncation=H)
        v.Pareto = c(v.Pareto,proc.time() -ptm)
        i = i +1
    }
    printf <- function(...) invisible(print(sprintf(...)))
    printf("ITM mean: %f", mean(v.IT))
    printf("ITM standard deviation: %f", sd(v.IT))
    printf("ARM mean: %f", mean(v.AR))
    printf("ARM standard deviation: %f", sd(v.AR))
    printf("built-in mean: %f", mean(v.Pareto))
    printf("built-in standard deviation: %f", sd(v.Pareto))
    
}

#compareTimes(65000,2000,0.25,2,3) # comentário para apagar: rodar quando entregar


## 1.2

## (a)

simexp=function(n,lam){
    x=0
    t=0
    while(t<n){
        t=t+1
        u=runif(1,0,1)
        x[t]=-log(u)/lam
    }
    x
}

sim.beta = function(n,a,b){
    if(a <= 0) stop("a must be a positive number")
    if(b <= 0) stop("b must be a positive number")
    if(n <= 0) stop("n must be a positive number")

    t=0
    Ya= 0 
    Yab= 0
    for(i in 1:(a+b)){
        Y = simexp(n,1)
        if(i <= a) {Ya = Ya + Y}
        Yab = Yab + Y
    }
    Ya/Yab
}

## (b)
set.seed(777)
a = 3
b = 1
m = 11000
simulation = sim.beta(m,a,b)
curveFun = function(x){(x^(a-1) * (1-x)^(b-1))/beta(a,b)}
p4 <- ggplot(data.frame(sample = simulation ), aes(x = sample)) +
     geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.1,
                   breaks=seq(0,1,0.1)) +
    geom_function(fun=curveFun)

plot(p4)

## (c)
set.seed(777)
sim.compare = rbeta(m,a,b)
quantile(simulation,type=1)
quantile(simulation,type=2)
quantile(simulation,type=3)
quantile(sim.compare,type=1)
quantile(sim.compare,type=2)
quantile(sim.compare,type=3)

## 2.1

## (a)

set.seed(1111)
f <- function(x){exp(-x^2/2)/sqrt(2*pi)}
i <- integrate(f, 2, Inf)$value; i

#b)
# Integral's integrand function : f(x) = exp(-x^2/2)/sqrt(2*pi) which is similar to the normal distribution.
# P(X > 2)
p <- 1-pnorm(2,0,1); p

## (c)

set.seed(1111)

g <- function(y){1/y^2 * exp(-(1/y)^2/2)/sqrt(2*pi)}

m=20000
x <- runif(m,0,1/2)

I.mc1 <- (1/2-0) * mean(g(x)); I.mc1
# 0.02301548
V.mc1 <- (1/2-0)^2 * var(g(x))/m; V.mc1
# 5.579652e-08

## (d)
set.seed(1111)

# Naive MC Estimator - Interval 0 to 1.
set.seed(1111)

g <- function(y){((exp(-1/(2*y^2)) * (1/y^2))- (exp((-(-(x+1))^2)/2)))/sqrt(2*pi)}
##g <- function(y){(2*exp((-(2/x)^2)/2))/(sqrt(2*pi) * ((x)^2))} não funciona bem na tabela
m=20000
x <- runif(m,0,1)

I.mc2 <- (1-0) * mean(g(x)); I.mc2
# 0.02408417
V.mc2 <- var(g(x))/m; V.mc2
# 1.543204e-06

## (e)

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

## (f)

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

##(g)
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



