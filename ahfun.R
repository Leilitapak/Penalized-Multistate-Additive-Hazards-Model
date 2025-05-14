get.coef <- function(cft, del, z) {
	ans <- .Call("get_coef", cft, as.integer(del), z, nrow(z))
	names(ans) <- c("b", "x", "v")
	ans
}

cd <- function(cft, del, z, method, lam, a, nlam=50, upto=0.4*nrow(z)) {
	ans <- get.coef(cft, del, z)
	b <- ans$b; x <- ans$x; v <- ans$v
	pen <- as.integer(switch(method, Lasso=0, SCAD=1, MCP=2, SICA=3, Enet=4))
	if (missing(lam)) {		# generate grid points for lambda
		lam.max <- max(abs(b)/v)
		lam <- lam.max*exp(seq(0, -log(nlam), length=nlam))
	}
	if (missing(a))			# default values for shape parameter
		a <- switch(method, Lasso=Inf, SCAD=3.7, MCP=3.7, SICA=c(1, 1e-4), Enet=0.05*(20:1))
	ans <- .Call("coord_descent", b, x, z, v, pen, lam, a, upto)
	sol <- ans[[1]]; nlam <- ans[[2]]
	list(sol=sol[, 1:nlam, ], lam=lam[1:nlam])
}

cv <- function(cft, del, z, method, lam, a, ind, nfold=10) {
	pen <- as.integer(switch(method, Lasso=0, SCAD=1, MCP=2, SICA=3, Enet=4))
	if (missing(a))
		a <- switch(method, Lasso=Inf, SCAD=3.7, MCP=3.7, SICA=c(1, 1e-4), Enet=0.05*(20:1))
	if (missing(ind)) {		# generate partition indices
		n <- length(cft); part <- integer(n)
		n1 <- sum(del); n0 <- n - n1
		part[del == 1] <- sample(rep(1:nfold, length=n1), n1)
		part[del == 0] <- sample(rep(1:nfold, length=n0), n0)
		ind <- split(1:n, part)
	}
	ans <- .Call("cross_valid", cft, as.integer(del), z, pen, lam, a, ind)
	names(ans) <- c("cv.err", "cv.se", "i", "j")
	ans
}
