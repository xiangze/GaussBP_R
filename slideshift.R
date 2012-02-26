expand <- function(img,r=1){
  row <- nrow(img)
  col <- ncol(img)
  nimg <- rbind(matrix(0,r,col+(2*r)),
                cbind(matrix(0,row,r),img,matrix(0,row,r)),
                matrix(0,r,col+(2*r)))
}

slidex <- function(inimg,d=1){
  row <- nrow(inimg)
  col <- ncol(inimg)
  d1 <- d+1
  img <-cbind(inimg[,d1:col],matrix(0,row,d))
}

slidexm <- function(inimg,d=1){
  row <- nrow(inimg)
  col <- ncol(inimg)
  d1 <- d+1
  img <-cbind(matrix(0,row,d),inimg[,1:(col-d)])
}

slidey <- function(inimg,d=1){
  row <- nrow(inimg)
  col <- ncol(inimg)
  d1 <- d+1
  img <- rbind(inimg[d1:row,],matrix(0,d,col))
}

slideym <- function(inimg,d=1){
  row <- nrow(inimg)
  col <- ncol(inimg)
  d1 <- d+1
  img <-rbind(matrix(0,d,col),inimg[1:(row-d),])
}


slideset4 <- function(inimg,d=1){
  row <- nrow(inimg)
  col <- ncol(inimg)
  xp <- slidex(inimg)
  xm <- slidexm(inimg)
  yp <- slidey(inimg)
  ym <- slideym(inimg)
  i4 <- list(xp=xp,xm=xm,yp=yp,ym=ym,c=inimg)
}

slideset4.4 <- function(inimg,d=1){
  hp <- slideset4(inimg$hp)
  hn <- slideset4(inimg$hn)
  vp <- slideset4(inimg$vp)
  vn <- slideset4(inimg$vn)
  i44 <-list(hp=hp,  hn=hn,  vp=vp,  vn=vn)
}
