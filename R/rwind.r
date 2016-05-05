# Wind Manipulation
# 2016-05-05 Ismail SEZEN
# sezenismail@gmail.com

#' rwind: A package for processing wind and related data.
#'
#' The rwind package provides four categories of important functions:
#' \code{\link{deg2comp}}, \code{\link{hd2uv}}, \code{\link{uv2hd}}
#' and \code{\link{rotax}}.
#'
#' @section rwind functions:
#' The rwind functions ...
#'
#' @docType package
#' @name rwind
NULL


index_array <- function(x, dim, value, drop = T) {
  indices <- rep(list(bquote()), length(dim(x)))
  indices[[dim]] <- value
  call <- as.call(c(
    list(as.name("["), quote(x)),
    indices,
    list(drop = drop)))
  eval(call)
}


windhelper_get_args <- function(...) {
  l <- list(...); x <- NULL
  if (length(l) == 1) {
    if (is.array(l[[1]])) x <- l[[1]]
    if (is.data.frame(l[[1]])) x <- as.matrix(l[[1]])
  }
  if (is.null(x)) x <- do.call(c, l)
  if (is.vector(x) && (length(x) %% 2) == 0)
    x <- array(x, dim = c(length(x) / 2, 2))
  return(x)
}


#' Wind direction to compass direction
#'
#' \code{deg2comp} converts horizontal wind speed and direction to uv
#' components of wind.
#'
#' @param x A numeric vector of wind directions in degrees.
#' @return Wind directions converted to compass directions
#'
#' @author Ismail SEZEN, \email{sezenismail@@gmail.com}
#'
#' @examples
#' deg2comp(runif(100, 0, 360))
#' x <- array(runif(1000, 0, 360), dim = c(10, 10, 10))
#' deg2comp(x)
#'
#' @export
deg2comp <- function(x, bins = c("N", "NNE", "NE", "ENE",
                                 "E", "ESE", "SE", "SSE",
                                 "S", "SSW", "SW", "WSW",
                                 "W", "WNW", "NW", "NNW")) {
  x[x == 0] <- 360
  x <- round((x * length(bins) / 360) + .5)
  x[1:length(x)] <- bins[x]
  return(x)
}


#' uv to horizontal speed & direction
#'
#' Convert uv components of wind to horizontal wind speed and direction.
#'
#' @param ... The names of the objects to be converted
#' @param deg Should the resulting wind directions are degree or radian?
#' @param reverse should be added 180 degree (or pi) to the wind directions?
#' @param drop Should be deleted the dimensions of an array which
#'        have only one level?
#' @return Speed and direction of the wind
#'
#' @author Ismail SEZEN, \email{sezenismail@@gmail.com}
#'
#' @examples
#' uv2hd(1, 1)
#' uv2hd(c(1, 1))
#' x <- array(runif(100, -5, 5), dim = c(100, 100, 2))
#' uv2hd(x)
#'
#' @export
uv2hd <- function(..., deg = T, reverse = T, drop = T) {
  x <- windhelper_get_args(...)
  comp_loc <- which(dim(x) == 2, arr.ind = T) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")
  dx <- dim(x); ldx <- length(dx)
  r <- aperm(x, c(comp_loc, (1:ldx)[-comp_loc])) # get uv dim to start

  num_pi <- if (deg) 180 else pi
  add_pi <- if (reverse) num_pi else 0
  u <- index_array(r, 1, 1)
  v <- index_array(r, 1, 2)
  h <- sqrt(colSums(r ^ 2))
  d <- ((atan2(u / h, v / h) * num_pi / pi) + add_pi)
  if (deg) d <- d  %% 360
  r <- array(c(h, d), dim = dx)

  r[is.nan(r)] <- 0
  dmn <- dimnames(x)
  if (is.null(dmn)) {
    dimnames(r)[[comp_loc]] <- c("h", "d")
  } else {
    dimnames(r) <- dmn
  }
  if (drop) r <- drop(r)
  return(r)
}


#' Horizontal wind speed & direction to uv
#'
#' Convert horizontal wind speed and direction to uv components of wind.
#'
#' @param ... The names of the objects to be converted
#' @param deg Should the resulting wind directions are degree or radian?
#' @param drop Should be deleted the dimensions of an array which
#'        have only one level?
#' @return uv components of the wind
#'
#' @author Ismail SEZEN, \email{sezenismail@@gmail.com}
#'
#' @examples
#' hd2uv(5, 30)
#' hd2uv(c(5, 30))
#' x <- matrix(c(runif(100, -5, 5), round(runif(100, 0, 360), 2)), ncol = 2,
#'             dimnames = list(NULL, c("h","d")))
#' hd2uv(x)
#'
#' @export
hd2uv <- function(..., deg = T, drop = T) {
  x <- windhelper_get_args(...)
  comp_loc <- which(dim(x) == 2, arr.ind = T) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents h and d")
  dx <- dim(x); ldx <- length(dx)
  r <- aperm(x, c(comp_loc, (1:ldx)[-comp_loc])) # get uv dim to start

  num_pi <- if (deg) 180 else pi
  h <- index_array(r, 1, 1)
  d <- index_array(r, 1, 2)
  d <- (d - num_pi) / num_pi
  u <- h * sinpi(d)
  v <- h * cospi(d)
  r <- array(c(u, v), dim = dx)

  r[is.nan(r)] <- 0
  dmn <- dimnames(x)
  if (is.null(dmn)) {
    dimnames(r)[[comp_loc]] <- c("u", "v")
  } else {
    dimnames(r) <- dmn
  }
  if (drop) r <- drop(r)
  return(r)
}


#' Rotate coordinate axis for uv wind components
#'
#' Rotates coordinate axes for uv components of the wind. uv values will be
#' re-calculated based on new axis. For instance, rotating axis 30 degrees
#' in clockwise means new and old positive y-directions will have 30 degrees
#'  between them.
#'
#'
#' @param ... The names of the objects to be converted.
#' @param alfa Rotation angle in degrees or radians.
#' @param deg TRUE if alfa is degree.
#' @param right TRUE if rotation is in clockwise.
#' @param drop Should be deleted the dimensions of an array which have
#'        only one level?
#' @return Rotated uv components
#'
#' @author Ismail SEZEN, \email{sezenismail@@gmail.com}
#'
#' @examples
#' rotax(1, 1)
#' rotax(1, 1, alfa = 120)
#' x <- array(runif(100, -5, 5), dim = c(100, 100, 2))
#' rotax(x)
#'
#' @export
rotax <- function(..., alfa = 45, deg = T, right = T, drop = T) {
  x <- windhelper_get_args(...)
  comp_loc <- which(dim(x) == 2, arr.ind = T) # find uv dim order
  if (length(comp_loc) == 0)
    stop("x does not have a dim = 2 represents u and v")

  u <- index_array(x, comp_loc, 1)
  v <- index_array(x, comp_loc, 2)
  a <- if (deg) alfa / 180 else alfa / pi
  c <- cospi(a); s <- sinpi(a)
  if (right) {
    un <- u * c - v * s
    vn <- u * s + v * c
  } else {
    un <- u * c + v * s
    vn <- v * c - u * s
  }
  dm <- dim(x); dn <- dimnames(x)
  x <- array(c(un, vn), dim = dm)
  dimnames(x) <- dn
  if (drop) x <- drop(x)
  return(x)
}
