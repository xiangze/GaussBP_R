library(RImageBook)

source("slideshift.R")

maxem <- 1000

imc <- readImage("sample.jpg")
ima <- imgRGB2Grey(imc)

n <-nrow(ima)
m <-ncol(ima)
l <- n*m

img <- expand(ima)
n1 <-nrow(img)
m1 <-ncol(img)

im <- slideset4(img)

mm <- matrix(0.001,n1,m1)
zm <- matrix(128.,n1,m1)

m4 <- list(xp=mm,xm=mm,yp=mm,ym=mm,c=mm)

mu <- list(hp=mm,hn=mm,vp=mm,vn=mm)
lambda <- list(hp=mm,hn=mm,vp=mm,vn=mm)

alpha.hz <- zm
sig.hz <- zm

dse.em <- 120
eps <- 100

#cat("initdone\n")

#cat("alpha.hz,sig.hz loop\n")

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

    #cat("t=",t,"\n")

    for (c in 1:nr){
      if(dse>eps){

        nit <- nit+1	

        muz <- slideset4.4(mu)
        lambdaz<- slideset4.4(lambda)

        #cat("inv \n")   
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

        #cat("ave.x \n")
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
        #cat("lambda \n")
        lambda$hp <- 1/(1/alpha+1/(isig2+lambdaz$hp$xm+lambdaz$vp$xm+lambdaz$vn$xm))
        lambda$hn <- 1/(1/alpha+1/(isig2+lambdaz$hn$xp+lambdaz$vp$xp+lambdaz$vn$xp))
        #cat("mu \n")
        
        mu$hp <- (im$xm*isig2+muz$hp$xm*lambdaz$hp$xm+muz$vn$xm*lambdaz$vn$xm
                    +muz$vp$xm*lambdaz$vp$xm)/
                      (isig2+lambdaz$hp$xm+lambdaz$vn$xm+lambdaz$vn$xm)

        mu$hn <- (im$xp*isig2+muz$hp$xp*lambdaz$hp$xp+muz$vn$xp*lambdaz$vn$xp
                    +muz$vp$xp*lambdaz$vp$xp)/
                      (isig2+lambdaz$hp$xm+lambdaz$vn$xm+lambdaz$vn$xm)

        #cat("inv \n")   
        ir11 <- isig2+alpha+lambdaz$hp$c+lambdaz$hn$c+lambdaz$vp$c
        ir22 <- isig2+alpha+lambdaz$hp$yp+lambdaz$hn$yp+lambdaz$vp$yp
        ir12 <- -alpha
        ir21 <- -alpha

        dd <- ir11*ir22-ir12*ir21
        rr11 <- ir22/dd;
        rr22 <- ir11/dd;
        rr12 <- -ir21/dd;
        rr21 <- ir12/dd;

        #cat("ave.y \n")
       
        ave.y1 <- (im$c*isig2+muz$hn$c*lambdaz$hn$c + muz$hp$c*lambdaz$hp$c
                   + muz$vp$c*lambdaz$vp$c+ muz$vn$c*lambdaz$vn$c)/
                     (isig2+lambdaz$hn$c+lambdaz$hp$c+lambdaz$vn$c+lambdaz$vp$c)

        ave.y2 <- (im$yp*isig2+muz$hn$yp*lambdaz$hn$yp + muz$hp$yp*lambdaz$hp$yp
                   + muz$vp$yp*lambdaz$vp$yp+ muz$vn$yp*lambdaz$vn$yp)/
                     (isig2+lambdaz$hn$yp+lambdaz$hp$yp+lambdaz$vn$yp+lambdaz$vp$yp)

        var.y1 <- r11+ave.y1*ave.y1
        var.y2 <- r22+ave.y2*ave.y2
        cov.y12 <- r12+ave.y1*ave.y2
        cor.y12 <- var.y1+var.y2-2*cov.y12

        lambda$vp <- 1/(1/alpha+1/(isig2+lambdaz$hp$ym+lambdaz$vp$ym+lambdaz$vn$ym))
        lambda$vn <- 1/(1/alpha+1/(isig2+lambdaz$hn$yp+lambdaz$vp$yp+lambdaz$vn$yp))

        mu$vp <- (im$ym*isig2+muz$hp$ym*lambdaz$hp$ym+muz$vn$ym*lambdaz$vn$ym
                    +muz$vp$ym*lambdaz$vp$ym)/
                      (isig2+lambdaz$hp$ym+lambdaz$vn$ym+lambdaz$vn$ym)

        mu$vn <- (im$yp*isig2+muz$hp$yp*lambdaz$hp$yp+muz$vn$yp*lambdaz$vn$yp
                  +muz$vp$yp*lambdaz$vp$yp)/
                    (isig2+lambdaz$hp$yp+lambdaz$vn$yp+lambdaz$vn$yp)

        dse.m <- apply(muz$hp$c-mu$hp,c(1,2),abs)
        +apply(muz$hn$c-mu$hn,c(1,2),abs)
        +apply(muz$vp$c-mu$vp,c(1,2),abs)
        +apply(muz$vn$c-mu$vn,c(1,2),abs)
        +apply(lambdaz$hp$c-lambda$hp,c(1,2),abs)
        +apply(lambdaz$hn$c-lambda$hn,c(1,2),abs)
        +apply(lambdaz$vp$c-lambda$vp,c(1,2),abs)
        +apply(lambdaz$vn$c-lambda$vn,c(1,2),abs)

        dse <- sum(dse.m)/l
#        muz <- mu
#        lambdaz <- lambda
#        fnit <- nit
#        fdse <- dse
      }#
    }#c

    cor.post.m <- cor.x12+cor.y12
    cor.post <- sum(cor.post.m)/(l*2)

    wk.m <-var.x1-2*im$c+ave.x1+im$c*im$c
    wk <- sum(wk.m)

    sig.hz <- sqrt(wk/l)
    alpha.hz <- 1/(2*cor.post)
    dse.em <- mean(abs(alpha-alpha.hz)+abs(1/(sig.h*sig.h)-1/(sig.hz*sig.hz)))
    cat("EM step:",t,"alpha_t=",alpha.hz,"sig_t=",sig.hz,"\n")
  }
}
zz <- ave.x1+0.5
writeImage(zz,"outimage.png")
cat("alpha_t=",alpha.hz,"sig_t=",sig.hz)
