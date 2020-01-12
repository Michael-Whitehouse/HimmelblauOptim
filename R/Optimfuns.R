
#' Himmelblau's function
#'
#' @param x coordinate
#' @param y coordinate
#'
#' @return Himmelblau's function at x,y
#' @export
#'
#' @examples
f<-function(x,y) {
  (x^2 + y -11)^2 + (x + y^2 - 7)^2
}





#' Himmelblau Gradient caluclator
#'
#' @param x x coordinate
#' @param y y coordinate
#'
#' @return value of the gradient at the given point
#' @export
#'
#' @examples
Himmelblau_grad <- function(x,y){
  c(4*x[1]*(x^2+y - 11) + 2*(x+y^2 - 7), 2*(x^2+y - 11) + 2*y*(x+y^2 - 7))
}




#' Gradient descent trajectory calculator
#'
#' @param x0      initial point
#' @param epsilon tolerance level
#' @param U       max step size
#'
#' @return an optimisation trajectory according to newton's method and given tolerance epsilon and max step size U
#' @export
#'
#' @examples
Himmelblau_grad_descent <- function(x0,epsilon,U){
  trajectory <- x0
  x <- x0
  gr <- Himmelblau_grad(x[1],x[2])
  while(norm(as.matrix(gr), type='F') > epsilon){

    t <- optim(0, function(t) f(x[1] -t*gr[1], x[2]-t*gr[2]),
               lower = 0, upper = U,  method =  'Brent' )$par
    x <- x - as.numeric(t)*gr
    gr <- Himmelblau_grad(x[1],x[2])
    trajectory <- cbind(trajectory, x)
  }
  return(trajectory)
}



#' Himmelblau Hessian calculator
#'
#' @param x x coordinate
#' @param y y coordinate
#'
#' @return value of the Hessian at the given point
#' @export
#'
#' @examples
Himmelblau_Hessian <- function(x,y){
  matrix(c(12*x^2+4*y -42, 4*(x+y),4*(x+y),6*y^2+2*x-12),2,2)
}


#' Newton's method trajectory calculator
#'
#' @param x0    initial point
#' @param epsilon  tolerance level
#'
#' @return an optimisation trajectory according to newton's method and given tolerance epsilon
#' @export
#'
#' @examples
Himmelblau_Newton <- function(x0,epsilon){
  trajectory <- x0
  x <- x0
  gr <- Himmelblau_grad(x[1],x[2])
  while(norm(as.matrix(gr), type='F') > epsilon){
    x <- x - solve(Himmelblau_Hessian(x[1],x[2]))%*%gr
    gr <- Himmelblau_grad(x[1],x[2])
    trajectory <- cbind(trajectory, x)
  }
  return(trajectory)
}

#' Title
#'
#' @param t trajectory
#'
#' @return a plot of the optimisation trajetory on the Himmelblau function xlim = ylim = c(-5,5)
#' @export
#'
#' @examples
plot_trajectory <- function(t){
  x1 <- seq(from = -5,to = 5, by =0.1) ## produce contour plot of Himmelblau's function
  contourmatrix <- outer(
    x1,
    x1,
    Vectorize(function(x,y) f(x,y)))
  contour(contourmatrix, nlevels = 50, drawlabels = F,xaxt='n',yaxt='n')



  points(t[1,1],t[2,1])
  for (i in c(1:(length(t[1,])-1))){
    segments(t[1,i],t[2,i],t[1,i+1],t[2,i+1], col = 'red')
  }
}
