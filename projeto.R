####################################################
##                        1.1                     ##
####################################################

pdf.pareto = function(x, alpha, L,H){
    (alpha * L^(alpha) * x^(-alpha-1))/(1 - (L/H)^(alpha))
}

## (b)
sim.IT = function(m,alpha, L, H){
    if(L>H) stop("Lower bound greater than Higher bound")
    else if(alpha <= 0) stop("Alpha must be a positive number")
    else if(L <= 0) stop("Lower Bound must be a positive number")
    else if(m <= 0) stop("n must be a positive number")
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



## (d) Since 15 samples is too small to make any conclusion, we cannot conclude that the code is not generating samples correctly. We can see that it has more samples in the beginning than in the end as we want so it can be likely that it is generating samples from the assumed distribution altough not conclusivly

library(ggplot2)
set.seed(2447)
alpha = 0.25
L = 2
H = 3
set.IT = sim.IT(12000,alpha,L,H)
par <- function(x){pdf.pareto(x,alpha,L,H)}
p <- ggplot(data.frame(it = set.IT[1:15]), aes(x = it)) +
     geom_histogram(aes(y = after_stat(density)),
                   binwidth = 0.2,
                   breaks=seq(L,H,0.1)) +
                   geom_function(fun=par)
plot(p)     

## (e)
sim.AR = function(n,alpha, L,H){
    if(L>H) stop("Lower bound greater than Higher bound")
    else if(alpha <= 0) stop("Alpha must be a positive number")
    else if(L <= 0) stop("Lower Bound must be a positive number")
    else if(n <= 0) stop("n must be a positive number")
    x = vector()
    x_rej = vector()
    y = vector()
    y_rej = vector()
    f =function (x){ (alpha * L^(alpha) * x^(-alpha-1))/(1 - (L/H)^(alpha)) }
    g = function(x){exp(-x)}
    h = function(x){f(x)/g(x)}
    M = optimize(h, c(L,H),maximum = T)$objective
    g.CDF=function(x){-exp(-x)+1}
    g.ITM = function(x){-log(1-x)}
    for(i in 1:n){
        u <- 1
        a <- 0
        while(u> a){
            x.c   <- g.ITM(runif(1,g.CDF(2),g.CDF(3))) # O g.CDF aqui serve para não estar a gerar valores desnecessários através da distribuição exponencial
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
    list(x = x,x_rej = x_rej,y = y, y_rej =y_rej, failRate= length(x_rej)/(length(x_rej) + length(x)))
}

## plot with rejected point and accepted points
alpha = 0.25
L = 2
H = 3
f =function (x){ (alpha * L^(alpha) * x^(-alpha-1))/(1 - (L/H)^(alpha)) }
g = function(x){exp(-x)}
h = function(x){f(x)/g(x)}
M = optimize(h, c(L,H),maximum = T)$objective
Mg= function(x){M*g(x)}
samples = sim.AR(12000,alpha,L,H)
plot(f,add,xlim=c(L,H),col="green")
curve(Mg,add=T,xlim=c(L,H))
points(samples$x,samples$y,pch=4,cex=1,lwd=1.5)
points(samples$x_rej,samples$y_rej,col=2,pch=4,cex=1,lwd=1.5)



## (g)
set.seed(2447)
n=12000
alpha=0.25
L=2
H=3
set.IT = sim.IT(n,alpha,L,H)
set.AR = sim.AR(n,alpha,L,H)$x
par(mfrow=c(1,2))
hist (set.IT,main="ITM Truncated Pareto",xlim=c(L,H),freq = F)
curve(pdf.pareto(x,alpha,L,H) ,add = T, xlim=c(L,H))
hist(set.AR,main="ARM Truncated Pareto", xlim=c(L,H),freq = F)
curve(pdf.pareto(x,alpha,L,H) ,add = T, xlim=c(L,H))    

## (f) d/dx h(x) = d/dx (0.25 * 2^(0.25) * x^(-1.25))/(1 - (2/3)^(0.25))/(e^-x) =
##     e^x (-3.85513/x^2.25 + 3.08411/x^1.25)
##     e^x (-3.85513/x^2.25 + 3.08411/x^1.25) = 0 => x ~= 1.25
##     Derivada positiva a partir de 1.25 logo está a crescer.
##     Queremos maior valor no intervalo [2,3]
##     Como a derivada só cresce o maior valor estará no 3
##     h(3) = 15.68958 <=> M = 15.68958
##     alpha = (1/M)h(x)

## (h)
compareTimes = function(m,nruns,alpha,L,H){
    library(Pareto)
    v.IT = vector()
    v.AR = vector()
    v.Pareto = vector()
    i=0
    while(i < nruns){
        print(i)
        ptm <- proc.time()
        z = sim.IT(m,alpha,L,H)
        v.IT = c(v.IT,proc.time() -ptm)
        ptm <- proc.time()
        z = sim.AR(m,alpha,L,H)
        v.AR = c(v.AR,proc.time() -ptm)
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

compareTimes(65000,2000,0.25,2,3)


####################################################
##                        1.2                     ##
####################################################

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


## (b) The simulation seems to follow the distribution correctly
set.seed(777)
a = 3
b = 1
m = 11000
simulation = sim.beta(m,a,b)
hist(simulation,freq=F)
curve((x^(a-1) * (1-x)^(b-1))/beta(a,b),add=T)

set.seed(777)
sim.compare = rbeta(m,a,b)
quantile(simulation,type=1)
quantile(simulation,type=2)
quantile(simulation,type=3)
quantile(sim.compare,type=1)
quantile(sim.compare,type=2)
quantile(sim.compare,type=3)
