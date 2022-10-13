plot.sshaped <- function(x, ...){
    #   some initial checkings
    if(!("sshaped" %in% class(x))) stop("Please supply an sshaped object.")
    
    if(x$shape=="sshaped") {shape <-"S-shaped fit"} else {shape <- "increasing convex fit"}
    plot(x$x,x$y,xlab="x",ylab="y", pch = 4, cex = 0.5, type = "p", main=shape)
    xorder = order(x$x)
    lines (x$x[xorder], x$fitted[xorder], ,lwd=3)
    if(x$shape=="sshaped") lines(c(x$inflection,x$inflection),c(min(x$y),max(x$y)),col="BLUE")
}
    

predict.sshaped <- function(object, xnew, ...){
    #   some initial checkings
    if(!("sshaped" %in% class(object))) stop("Please supply an sshaped object.")
    if((!is.vector(xnew))||(!is.numeric(xnew)))
        stop("Please check the type of the new preditor - need to be numeric vectors")
    if (missing(xnew)) return(object$fitted) 
    xorder = order(object$x)
    x = object$x[xorder]
    f = object$fitted[xorder]
    n = length(x)
    if (n==1) ret<-rep(f[1],length(xnew))
    else {
        predict_single <-function(singlex){
            if(singlex < x[1]) fvalue <-(singlex-x[1])*(f[2]-f[1])/(x[2]-x[1])+f[1]
            else if (singlex > x[n]) fvalue <-(singlex-x[n])*(f[n]-f[n-1])/(x[n]-x[n-1])+f[n]
            else fvalue <- approx(x,f,singlex)$y
        }
        ret<-sapply(xnew, predict_single)
    }
    return(ret)
}

