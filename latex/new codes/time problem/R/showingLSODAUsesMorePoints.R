
# this file contains the code for the tolerance study of lsoda for the 
# time problem with and without event detection

library(deSolve)
library(ggplot2)
library(plyr)

## defining the required functions

### no event detection
model_no_event <- function(t, state, parms) {
    S <- state[1]
    E <- state[2]
    I <- state[3]
    R <- state[4]

    N <- 37.741 * (10^6)
    alpha <- 1.0/8.0
    beta <- 0.9
    if (t > 27) {
        beta <- 0.005
    }
    gamma <- 0.06
    mu <- 0.01/365

    dSdt <- mu*N -mu*S - (beta/N)*I*S
    dEdt <- (beta/N)*I*S - alpha*E - mu*E
    dIdt <- alpha*E - gamma*I -mu*I
    dRdt <- gamma*I - mu*R
    
    list(c(dSdt, dEdt, dIdt, dRdt))
}


experiment_no_event <- function(by) {
    t0 <- 0
    tspan <- seq(from=0, to=95, by=by)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0 <- c(S0, E0, I0, R0)

    out <- ode(func=model_no_event, y=y0, times=tspan, parms=NULL, method="lsoda")
    diag <- diagnostics(out)

    colnames(out) <- c("time",
                        paste("S", by, sep="_"),
                        paste("E", by, sep="_"),
                        paste("I", by, sep="_"),
                        paste("R", by, sep="_")

    )
    out
}

### with event functions
model_before <- function(t, state, parms) {
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
    
    list(c(dSdt, dEdt, dIdt, dRdt))
}

model_after <- function(t, state, parms) {
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
    
    list(c(dSdt, dEdt, dIdt, dRdt))
}

experiment_with_event <- function(by) {
    t0 <- 0
    tspan_before <- seq(from=0, to=27, by=by)
    tspan_after <- seq(from=27, to=95, by=by)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0_before <- c(S0, E0, I0, R0)

    out_before <- ode(func=model_before, y=y0_before, times=tspan_before, 
                      parms=NULL, method="lsoda")
    y0_after <- (out_before[nrow(out_before),])[!(colnames(out_before) %in% c("time"))]
    out_after <- ode(func=model_after, y=y0_after, times=tspan_after, 
                         parms=NULL, method="lsoda")
    out <- rbind(out_before, out_after)
    colnames(out) <- c("time",
                        paste("S", by, sep="_"),
                        paste("E", by, sep="_"),
                        paste("I", by, sep="_"),
                        paste("R", by, sep="_")
    )
    out
}


lsoda_no_event_1 <- experiment_no_event(1)
lsoda_no_event_3 <- experiment_no_event(3)
lsoda_no_event_5 <- experiment_no_event(5)

plot(lsoda_no_event_1[, 1], lsoda_no_event_1[, 4], type="l", col="black", lwd=1.2)
lines(lsoda_no_event_3[, 1], lsoda_no_event_3[, 4], col="blue", lwd=1.2)
lines(lsoda_no_event_5[, 1], lsoda_no_event_5[, 4], col="red", lwd=1.2)
legend(1, 25000, 
       legend=c("lsoda_no_event_1", "lsoda_no_event_3", "lsoda_no_event_5"),
       col=   c("black", "blue", "red"), lty=1:2, cex=0.8,
       title="Line types", text.font=4, bg='lightblue')

# =================================================

#lsoda_with_event_1 <- experiment_with_event(1)
#lsoda_with_event_3 <- experiment_with_event(3)
#lsoda_with_event_5 <- experiment_with_event(5)

#plot(lsoda_with_event_1[, 1], lsoda_with_event_1[, 4], type="l", col="black", lwd=1.2)
#lines(lsoda_with_event_3[, 1], lsoda_with_event_3[, 4], col="blue", lwd=1.2)
#lines(lsoda_with_event_5[, 1], lsoda_with_event_5[, 4], col="red", lwd=1.2)
#legend(1, 25000, 
#       legend=c("lsoda_with_event_1", "lsoda_with_event_3", "lsoda_with_event_5"),
#       col=   c("black", "blue", "red"), lty=1:2, cex=0.8,
#       title="Line types", text.font=4, bg='lightblue')
