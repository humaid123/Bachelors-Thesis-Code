
library(deSolve)
library(ggplot2)
library("plyr")
library(hash)
library(jsonlite)

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

experiment <- function(method, atol, rtol) {
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
                       rootfun=root_10000, atol=atol, rtol=rtol)
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
    list(res=res, nfev=nfev, time=as.numeric(time))
}

ans <- hash()
tolerances <- c(1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-15)
for (tolerance in tolerances) {
    methods <- c("lsoda", "radau", "bdf", "adams")
    for (method in methods) {
         ans[[paste(method, round(log10(tolerance)), sep="_")]] <- experiment(method, tolerance, tolerance)
    }
}
ans <- as.list.hash(ans)
write_json(ans, "tolerance_event_R.json")
