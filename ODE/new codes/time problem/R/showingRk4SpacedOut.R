# R's rk4 is the classical butcher: see call_rk4.c main loop =>
# it evaluates the function at t, t+0.5h, t+0.5h and t+h

# we can also look at call_euler.c => the code does not estimate the error
# it just brute-force from left to right

# the timestep it uses is the spacings in times:
# tt is the times we want to sample
# time_step = dt = tt[i+1] - tt[i]

# So if we sample it less frequently, it should mess up

library(deSolve)

model <- function(t, state, parms) {
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

t0 <- 0
tspan <- seq(from=0, to=95, by=1)
N <- 37.741 * (10^6)
E0 <- 103
I0 <- 1
R0 <- 0
S0 <- N - (E0 + I0 + R0)

y0 <- c(S0, E0, I0, R0)

# lsoda is the default function which we know gives the ANS
lsoda <- ode(func=model, y=y0, times=tspan, parms=NULL)

# tring out rk4 with a worse tspan
tspan <- seq(from=0, to=95, by=1)
rk4_tstep_1 <- ode(func=model, y=y0, times=tspan,  parms=NULL, method="rk4")

tspan <- seq(from=0, to=95, by=2)
rk4_tstep_2 <- ode(func=model, y=y0, times=tspan,  parms=NULL, method="rk4")

tspan <- seq(from=0, to=95, by=5)
rk4_tstep_5 <- ode(func=model, y=y0, times=tspan,  parms=NULL, method="rk4")

tspan <- seq(from=0, to=95, by=7)
rk4_tstep_7 <- ode(func=model, y=y0, times=tspan,  parms=NULL, method="rk4")

plot(rk4_tstep_2[, 1], rk4_tstep_2[, 4], type="l", col="black", lwd=2,
    xlab="time", 
    ylab="I(t)", 
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
    cex.sub = 1.5
)
lines(rk4_tstep_5[, 1], rk4_tstep_5[, 4], col="orange", lwd=2)
lines(rk4_tstep_7[, 1], rk4_tstep_7[, 4], col="green", lwd=2)
lines(rk4_tstep_1[, 1], rk4_tstep_1[, 4], col="blue", lwd=2)
lines(lsoda[, 1], lsoda[, 4], type="l", col="red", lwd=2)
legend(1, 30000, 
       legend=c("space=2", "space=5", "space=7", "space=1", "lsoda"),
       col=   c("black", "orange", "green", "blue", "red"), lty=1, lwd=2,
       title="legend", 
       cex=1.1, text.font=4)



### EVENT WITH RK4
library(deSolve)

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

experiment <- function(method, step) {
    t0 <- 0
    tspan_before <- c(seq(from=0, to=27, by=step), 27)
    tspan_after <- seq(from=27, to=95, by=step)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0_before <- c(S0, E0, I0, R0)

    out_before <- ode(func=model_before, y=y0_before, times=tspan_before, parms=NULL, method=method)

    y0_after <- (out_before[nrow(out_before),])[!(colnames(out_before) %in% c("time"))]
    out_after <- ode(func=model_after, y=y0_after, times=tspan_after, parms=NULL, method=method)

    out <- rbind(out_before, out_after)
    colnames(out) <- c("time",
                        paste("S", method, step, sep="_"),
                        paste("E", method, step, sep="_"),
                        paste("I", method, step, sep="_"),
                        paste("R", method, step, sep="_")
    )
    out
}

# lsoda is the default function which we know gives the ANS
lsoda <- experiment("lsoda", 1)

# trying out rk4 with a worse tspan
rk4_tstep_1 <- experiment("rk4", 1)

rk4_tstep_2 <- experiment("rk4", 2)

rk4_tstep_5 <- experiment("rk4", 5)

rk4_tstep_7 <- experiment("rk4", 7)

plot(rk4_tstep_2[, 1], rk4_tstep_2[, 4], type="l", col="black", lwd=2,
    xlab="time", 
    ylab="I(t)", 
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
    cex.sub = 1.5
)
lines(rk4_tstep_5[, 1], rk4_tstep_5[, 4], col="orange", lwd=2)
lines(rk4_tstep_7[, 1], rk4_tstep_7[, 4], col="green", lwd=2)
lines(rk4_tstep_1[, 1], rk4_tstep_1[, 4], col="blue", lwd=2)
lines(lsoda[, 1], lsoda[, 4], type="l", col="red", lwd=2)
legend(1, 25000, 
       legend=c("space=2", "space=5", "space=7", "space=1", "lsoda"),
       col=   c("black", "orange", "green", "blue", "red"), lty=1, lwd=2,
       title="legend", 
       cex=1.1, text.font=4)
