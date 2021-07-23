library(deSolve)
library(ggplot2)
library("plyr")

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
    time <- 0
    nfev <- 0

    res <- matrix(nrow=0, ncol=5)
    t_initial <- t0
    y_initial <- y0
    measures_implemented <- FALSE


    while (t_initial < 180) {
        tspan = seq(from=t_initial, to=180, by=1)
        out <- 1
        if (measures_implemented) {
            tic <- Sys.time()   
            out <- ode(func=model_with_measures, y=y_initial, 
                       times=tspan, parms=NULL, method=method, 
                       rootfun=root_10000)
            toc <- Sys.time()   
            time <- time + (toc-tic)
            diag <- diagnostics(out)
            nfev <- nfev + diag$istate[3]
            measures_implemented <- FALSE
        } else {
            tic <- Sys.time()
            out <- ode(func=model_no_measures, y=y_initial, 
                       times=tspan, parms=NULL, method=method,
                       rootfun=root_25000)
            toc <- Sys.time()   
            time <- time + (toc-tic)
            diag <- diagnostics(out)
            nfev <- nfev + diag$istate[3]
            measures_implemented <- TRUE
        }
        
        t_change <- length(out[, 1])
        t_initial <- round(out[t_change, 1]) #round(out[t_change, 1])
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
    list(res=res, nfev=nfev, time=time)
}

out <- experiment("lsoda")
lsoda_res <- out$res
lsoda_nfev <- out$nfev
lsoda_time <- out$time
#daspk_res <- experiment("daspk") => no root finding
#euler_res <- experiment("euler") => no root finding
#rk4_res <- experiment("rk4")     => no root finding
#ode45_res <- experiment("ode45") => no root finding
out <- experiment("radau")
radau_res <- out$res
radau_nfev <- out$nfev
radau_time <- out$time

out <- experiment("bdf")
bdf_res <- out$res
bdf_nfev <- out$nfev
bdf_time <- out$time

out <- experiment("adams")
adams_res <- out$res
adams_nfev <- out$nfev
adams_time <- out$time

#p <- ggplot(as.data.frame(res), aes(x=time))
#p <- p + geom_line(aes(y=E_root), color="green")
#p <- p + geom_line(aes(y=I_root), color="blue")
#p

final <- as.data.frame(lsoda_res)
final <- join(final, as.data.frame(radau_res), by="time", type="left")
final <- join(final, as.data.frame(bdf_res), by="time", type="left")
final <- join(final, as.data.frame(adams_res), by="time", type="left")

colors <- c("lsoda"="blue", "radau"="orange", "bdf"="brown", "adams"="red")
p <- ggplot(as.data.frame(final), aes(x=time)) + 
    geom_line(aes(y=E_lsoda, color="lsoda"), size=1.2) +
    geom_line(aes(y=E_radau, color="radau"), size=1.2) +
    geom_line(aes(y=E_bdf, color="bdf"), size=1.2)     + 
    geom_line(aes(y=E_adams, color="adams"), size=1.2) +
    labs(x="time", y="E", color="legend") +
    scale_color_manual(values=colors)
p

m <- rbind(
           c("lsoda", lsoda_nfev, lsoda_time), 
           c("radau", radau_nfev, radau_time), 
           c("bdf", bdf_nfev, bdf_time), 
           c("adams", adams_nfev, adams_time) 
           )
m

# lsoda
#h1 <- ggplot(as.data.frame(final), aes(x=time))
#h1 <- h1 + geom_line(aes(y=E_lsoda), color="green")
#h1 <- h1 + geom_line(aes(y=I_lsoda), color="blue")
#h1

# radau
#h2 <- ggplot(as.data.frame(final), aes(x=time))
#h2 <- h2 + geom_line(aes(y=E_radau), color="green")
#h2 <- h2 + geom_line(aes(y=I_radau), color="blue")
#h2

# bdf
#h3 <- ggplot(as.data.frame(final), aes(x=time))
#h3 <- h3 + geom_line(aes(y=E_bdf), color="green")
#h3 <- h3 + geom_line(aes(y=I_bdf), color="blue")
#h3

# adams
#h4 <- ggplot(as.data.frame(final), aes(x=time))
#h4 <- h4 + geom_line(aes(y=E_adams), color="green")
#h4 <- h4 + geom_line(aes(y=I_adams), color="blue")
#h4
