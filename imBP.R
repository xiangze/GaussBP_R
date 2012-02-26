library(RImageBook)

source("slideshift.R")

maxem <- 10

imc <- readImage("sample.jpg")
ima <- imgRGB2Grey(imc)

n <-nrow(ima)
m <- ncol(ima)

img <- expand(ima)
n1 <- n+1
m1 <- m+1

im <- slideset4(img)

mm <- matrix(0.001,n1,m1)
zm <- matrix(128.,n1,m1)

m4 <- list(xp=mm,xm=mm,yp=mm,ym=mm)

mu <- list(hp=m4,hn=m4,vp=m4,vn=m4)
lambda <- list(hp=m4,hn=m4,vp=m4,vn=m4)

#muz <- mu
#lambdaz <- mu

ave <- list(xp=zm,xm=zm,yp=zm,ym=zm)
var <- list(xp=zm,xm=zm,yp=zm,ym=zm)

alpha.hz <- zm
sig.hz <- zm

dse.em <- 120
eps <- 100

for (t in 1:maxem){
  if(dse.em> 0.00001){

    alpha.h <- alpha.hz
    sig.h <- sig.hz

    alpha <- alpha.h
    sig <- sig.h
    isig2 <- 1/(sig*sig)

    nit <- 0
    nr <- 500
    dse <- 10000000

    for (c in 1:nr){
      if(dse>eps){
        nit <- nit+1	

        muz <- slideset4.4(mu)
        lambdaz<- slideset4.4(lambda)

        cat("inv \n")   
        ir11 <- isig2+alpha+lambdaz$hp$c+lambdaz$vp$c+lambdaz$vn$c
        ir22 <- isig2+alpha+lambdaz$hp$xp+lambdaz$vp$xp+lambdaz$vn$xp
        ir12 <- -alpha
        ir21 <- -alpha 
#inv
        dd <- ir11*ir22-ir12*ir21
        r11 <- ir22/dd
        r22 <- ir11/dd
        r12 <- -ir21/dd
        r21 <- ir12/dd

        cat("ave \n")
        ave.x1 <- (im$c*isig2+muz$hn$c*lambdaz$hn$c + muz$hp$c*lambdaz$hp$c
                   + muz$vp$c*lambdaz$vp$c+ muz$vn$c*lambdaz$vn$c)/
                     (isig2+lambdaz$hn$c + lambdaz$hp$c +lambdaz$vp$c + lambdaz$vn$c)

        ave.x2 <- (im$xp*isig2+muz$hn$xp*lambdaz$hn$xp + muz$hp$xp*lambdaz$hp$xp
                   + muz$vp$xp*lambdaz$vp$xp+ muz$vn$xp*lambdaz$vn$xp)/
                     (isig2+lambdaz$hn$xp + lambdaz$hp$xp +lambdaz$vp$xp+ lambdaz$vn$xp)

        var.x1 <- r11+ave.x1*ave.x1
        var.x2 <- r22+ave.x2*ave.x2
        cov.x12 <- r12+ave.x1*ave.x2
        cor.x12 <- var.x1+var.x2-2*cov.x12

        lambda$hp <- 1/(1/alpha+1/(isig2+labmdaz$hp$xm+labmdaz$vp$xm+labmdaz$vm$xm))
        lambda$hn <- 1/(1/alpha+1/(isig2+labmdaz$hn$xp+labmdaz$vp$xp+labmdaz$vm$xp))

        mu$hp <- (im$xm*isig2+muz$hp$xm*lambdaz$hp$xm+muz$vn$xm*lambdaz$vn$xm
                    +muz$vp$xm*lambdaz$vp$xm)/
                      (isig2+lambdaz$hp$xm+lambdaz$vn$xm+lambdaz$vn$xm)

        mu$hn <- (im$xp*isig2+muz$hp$xp*lambdaz$hp$xp+muz$vn$xp*lambdaz$vn$xp
                    +muz$vp$xp*lambdaz$vp$xp)/
                      (isig2+lambdaz$hp$xm+lambdaz$vn$xm+lambdaz$vn$xm)

        ir11 <- isig+alpha+lambdaz$hp$c+lambdaz$hn$c+lambdaz$vp$c
        ir22 <- isig+alpha+lambdaz$hp$yp+lambdaz$hn$yp+lambdaz$vp$yp
        ir12 <- -alpha
        ir21 <- -alpha
#inv
        dd <- ir11*ir22-ir12*ir21
        rr11 <- ir22/dd;
        rr22 <- ir11/dd;
        rr12 <- -ir21/dd;
        rr21 <- ir12/dd;

        ave.y1 <- (im$c*isig2+muz$hn$c*lambdaz$hn$c + muz$hp$c*lambdaz$hp$c
                   + muz$vp$c*lambdaz$vp$c+ muz$vn$c*lambdaz$vn$c)/
                     (isig2+lambdaz$hn$c+lambdaz$hp$c+lambdaz$vn$c+lambdaz$vp$c)

        ave.y2 <- (im$yp*isig2+muz$hn$yp*lambdaz$hn$yp + muz$hp$yp*lambdaz$hp$yp
                   + muz$vp$yp*lambdaz$vp$yp+ muz$vn$yp*lambdaz$vn$yp)/
                     (isig2+lambdaz$hn$yp+lambdaz$hp$yp+lambdaz$vn$yp+lambdaz$vp$yp)

        var.y1 <- r11+ave.y1*ave.y1
        var.y2 <- r22+ave.y2*ave.y2
        cov.y12 <- r12+ave.y1*ave.xy
        cor.y12 <- var.y1+var.y2-2*cov.y12

        lambda$vp <- 1/(1/alpha+1/(isig2+labmdaz$hp$ym+labmdaz$vp$ym+labmdaz$vm$ym))
        lambda$vn <- 1/(1/alpha+1/(isig2+labmdaz$hn$yp+labmdaz$vp$yp+labmdaz$vm$yp))

        mu.x$vp <- (im.ym*isig2+muz.ym$hp*lambdaz.ym$hp+muz.ym$vn*lambdaz.ym$vn+muz.ym$vp*lambdaz.ym$vp)/(isig2+lambdaz.ym$hp+lambdaz.ym$vn+lambdaz.ym$vn)

        mu.x$vn <- (im.y*isig2+muz.y$hp*lambdaz.y$hp+muz.y$vn*lambdaz.y$vn+muz.y$vp*lambdaz.y$vp)/(isig2+lambdaz.y$hp+lambdaz.y$vn+lambdaz.y$vn)

        dse.m <- apply(muz$hp$c-mu$hp$c,c(1,2),abs)
        +apply(muz$hn$c-mu$hn$c,c(1,2),abs)
        +apply(muz$vp$c-mu$vp$c,c(1,2),abs)
        +apply(muz$vn$c-mu$vn$c,c(1,2),abs)
        +apply(lambdaz$hp$c-lambda$hp$c,c(1,2),abs)
        +apply(lambdaz$hn$c-lambda$hn$c,c(1,2),abs)
        +apply(lambdaz$vp$c-lambda$vp$c,c(1,2),abs)
        +apply(lambdaz$vn$c-lambda$vn$c,c(1,2),abs)

        dse <- sum(dse.m)/l

#        muz <- mu
#        lambdaz <- lambda
        cor.post.m <- cor.x12+cor.x12
        cor.post <- sum(cor.post.m)/(l*l)

        zz <- ave.x1+0.5

        fnit <- nit
        fdse <- dse
      }#
    }#c

    wk.m <-var.x1-2*im+ave.x1+im*im
    wk <- sum(wk.m)

    sigmahz <- sqrt(wk/l)
    alpha.hz <- 1/(2*cor.post)
    dse.em <- abs(alpha-alpha.hz)+abs(1/(sig.h*sig.h)-1/(sig.hz*sig.hz))
    print(c("EM step:",t,"alpha_t=",alpha.hz,"sig_t=",sig.hz))

  }
}

imwrite(zz,"outimage.png")
print(c("alpha_t=",alpha.hz,"sig_t=",sig.hz))
