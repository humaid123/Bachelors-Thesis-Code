

library(deSolve)
library(ggplot2)
library(plyr)
library(hash)
library(jsonlite)

measures_implemented <- 0
direction <- "up"

model_with_if <- function(t, state, parms) {
    S <- state[1]
    E <- state[2]
    I <- state[3]
    R <- state[4]

    N <- 37.741 * (10^6)
    alpha = 1.0/8.0
    if (.GlobalEnv$direction == "up") {
        if (E > 25000) {
            .GlobalEnv$measures_implemented <- 1
            .GlobalEnv$direction <- "down"
        }
    } else {
        if (E < 10000) {
            .GlobalEnv$measures_implemented <- 0
            .GlobalEnv$direction <- "up"
        }
    }

    beta <- 0.9
    if (.GlobalEnv$measures_implemented == 1) {
        beta <- 0.005
    }
    gamma = 0.06
    mu = 0.01/365

    dSdt <- mu*N - mu*S - (beta/N)*I*S
    dEdt <- (beta/N)*I*S - alpha*E - mu*E
    dIdt <- alpha*E - gamma*I - mu*I
    dRdt <- gamma*I - mu*R

    list(c(dSdt, dEdt, dIdt, dRdt))
}

experiment_with_if <- function(method, atol, rtol) {
    t0 <- 0
    tspan <- seq(from=0, to=180, by=1)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0 <- c(S0, E0, I0, R0)


    .GlobalEnv$measures_implemented <- 0
    .GlobalEnv$direction <- "up"
    # lsoda is the default function
    tic <- Sys.time()
    out <- ode(func=model_with_if, y=y0, times=tspan, parms=NULL, 
               method=method, atol=atol, rtol=rtol)
    toc <- Sys.time()
    diag <- diagnostics(out)
    colnames(out) <- c("time",
                        paste("S", method, sep="_"),
                        paste("E", method, sep="_"),
                        paste("I", method, sep="_"),
                        paste("R", method, sep="_")
    )
    list(res=out, nfev=diag$istate[3], time=as.numeric(toc-tic))
}

ans <- hash()
tolerances <- c(1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-15)
for (tolerance in tolerances) {
    methods <- c("lsoda", "daspk", "euler", "rk4", "ode45", 
                 "radau", "bdf", "adams")
    for (method in methods) {
         ans[[paste(method, round(log10(tolerance)), sep="_")]] <- experiment_with_if(method, tolerance, tolerance)
    }
}
ans <- as.list.hash(ans)
write_json(ans, "trial.json")
