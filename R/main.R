cvxreg <- function(x, y) {
    #   some initial checkings
    if((!is.vector(x))||(!is.vector(y))||(!is.numeric(x))||(!is.vector(y)))
        stop("Please check the type of your inputs - they need to be numeric vectors")
    if(length(x)!=length(y))
        stop("Please check your inputs - they need to be of same length")
    if(length(unique(x))!=length(x))
        stop("Please check the design points - in this version, only distinct covariates are supported")
    
    #   call function coded via rcpp
    output <- .Call(`_Sshaped_cvxreg`, x, y)

    #   produce output in a sshaped class object
    ret <- list()	  
    ret$x <- x
    ret$y <- y
    ret$fitted <- output[[3]][order(order(x))]
    ret$rss <- sum((ret$fitted-ret$y)^2)
    ret$inflection <- max(x) 
    ret$shape <- "convex"
    class(ret) <- "sshaped"
    return(ret)
}

sshapedreg <- function(x, y) {
    #   some initial checkings
    if((!is.vector(x))||(!is.vector(y))||(!is.numeric(x))||(!is.vector(y)))
        stop("Please check the type of your inputs - they need to be numeric vectors")
    if(length(x)!=length(y))
        stop("Please check your inputs - they need to be of same length")
    if(length(unique(x))!=length(x))
        stop("Please check the design points - in this version, only distinct covariates are supported")
    
    #   call function coded via rcpp
    output<- .Call(`_Sshaped_sshapedreg`, x, y)
    
    
    #   produce output in a sshaped class object
    ret <- list()	  
    ret$x <- x
    ret$y <- y
    ret$fitted <- output[[2]][order(order(x))]
    ret$rss <- sum((ret$fitted-ret$y)^2)
    ret$inflection <-  output[[1]]
    ret$shape <- "sshaped"
    class(ret) <- "sshaped"
    return(ret)
}
