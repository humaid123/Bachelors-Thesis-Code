library(deSolve)
library(ggplot2)
library(plyr)

model_no_measures <- function(t, state, parms) {
    S <- state[1]
    E <- state[2]
    I <- state[3]
    R <- state[4]

    N <- 37.741 * (10^6)
    alpha <- 1.0/8.0
    beta <- 0.9
    gamma <- 0.06
    mu <- 0.01/365

    dSdt <- mu*N -mu*S - (beta/N)*I*S
    dEdt <- (beta/N)*I*S - alpha*E - mu*E
    dIdt <- alpha*E - gamma*I -mu*I
    dRdt <- gamma*I - mu*R
    res <- list(c(dSdt, dEdt, dIdt, dRdt))
    res
}

root_25000 <- function(t, state, parms) {
    E <- state[2]
    E - 25000
}

model_with_measures <- function(t, state, parms) {
    S <- state[1]
    E <- state[2]
    I <- state[3]
    R <- state[4]

    N <- 37.741 * (10^6)
    alpha <- 1.0/8.0
    beta <- 0.005
    gamma <- 0.06
    mu <- 0.01/365

    dSdt <- mu*N -mu*S - (beta/N)*I*S
    dEdt <- (beta/N)*I*S - alpha*E - mu*E
    dIdt <- alpha*E - gamma*I -mu*I
    dRdt <- gamma*I - mu*R
    res <- list(c(dSdt, dEdt, dIdt, dRdt))
    res
}

root_10000 <- function(t, state, parms) {
    E <- state[2]
    E - 10000
}

experiment <- function(method) {
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0 <- c(S0, E0, I0, R0)
    t0 <- 0
    nfev <- 0

    res <- matrix(nrow=0, ncol=5)
    t_initial <- t0
    y_initial <- y0
    measures_implemented <- FALSE

    # need to use 179 for t_span to have at least one value
    while (t_initial < 180) {
        tspan = seq(from=t_initial, to=180, by=1)
        if (length(tspan) == 1) {
            break
        }
        out <- 1
        if (measures_implemented) {
            out <- ode(func=model_with_measures, y=y_initial, 
                       times=tspan, parms=NULL, method=method, 
                       rootfun=root_10000)
            diag <- diagnostics(out)
            nfev <- nfev + diag$istate[3]
            measures_implemented <- FALSE
        } else {
            out <- ode(func=model_no_measures, y=y_initial, 
                       times=tspan, parms=NULL, method=method,
                       rootfun=root_25000)
            diag <- diagnostics(out)
            nfev <- nfev + diag$istate[3]
            measures_implemented <- TRUE
        }
        
        t_change <- length(out[, 1])
        t_initial <- (out[t_change, 1]) #round(out[t_change, 1])
        y_initial <- out[t_change, ]
        y_initial <- y_initial[!(names(y_initial) %in% c("time"))]
        res <- rbind(res, out)
    }

    colnames(res) <- c("time", 
                        paste("S", method, sep="_"),
                        paste("E", method, sep="_"),
                        paste("I", method, sep="_"),
                        paste("R", method, sep="_")
                      )
    list(res=res, nfev=nfev)
}

out <- experiment("lsoda")
lsoda_res <- out$res
lsoda_nfev <- out$nfev
#daspk_res <- experiment("daspk") => no root finding
#euler_res <- experiment("euler") => no root finding
#rk4_res <- experiment("rk4")     => no root finding
#ode45_res <- experiment("ode45") => no root finding
out <- experiment("radau")
radau_res <- out$res
radau_nfev <- out$nfev

out <- experiment("bdf")
bdf_res <- out$res
bdf_nfev <- out$nfev

out <- experiment("adams")
adams_res <- out$res
adams_nfev <- out$nfev

m <- rbind(
    c("lsoda", lsoda_nfev), 
    c("radau", radau_nfev), 
    c("bdf", bdf_nfev), 
    c("adams", adams_nfev) 
)
m

plot(lsoda_res[, 1], lsoda_res[, 3], type="l", col="black", lwd=2, 
    xlab="time", 
    ylab="E(t)", 
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
    cex.sub = 1.5
)
lines(radau_res[, 1], radau_res[, 3], col="orange", lwd=2)
lines(adams_res[, 1], adams_res[, 3], col="green", lwd=2)
lines(bdf_res[, 1], bdf_res[, 3], col="blue", lwd=2)
legend(-6, 25000, 
       legend=c("lsoda", "radau", "adams", "bdf"),
       col=   c("black", "orange", "green", "blue"), lty=1, lwd=2,
       title="Legend", text.font=4, cex=1.1)
