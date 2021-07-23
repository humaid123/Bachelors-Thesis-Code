
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

# ode45 spaced out
tspan <- seq(from=0, to=95, by=1)
tstep_1 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_1)
count <- c("tstep_1", as.integer(res$istate[3]))
m <- count


tspan <- seq(from=0, to=95, by=3)
tstep_3 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_3)
count <- c("tstep_3", as.integer(res$istate[3]))
m <- rbind(m, count)


tspan <- seq(from=0, to=95, by=5)
tstep_5 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_5)
count <- c("tstep_5", as.integer(res$istate[3]))
m <- rbind(m, count)

tspan <- seq(from=0, to=95, by=7)
tstep_7 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_7)
count <- c("tstep_7", as.integer(res$istate[3]))
m <- rbind(m, count)


plot(tstep_7[, 1], tstep_7[, 4], type="l", col="black", lwd=2)
lines(tstep_5[, 1], tstep_5[, 4], type="l", col="brown", lwd=2)
lines(tstep_3[, 1], tstep_3[, 4], col="orange", lwd=2)
lines(tstep_1[, 1], tstep_1[, 4], col="blue", lwd=2)
lines(lsoda[, 1], lsoda[, 4], type="l", col="red", lwd=2)
legend(1, 30000, 
       legend=c("space=7", "space=5", "space=3", "space=1", "ANS"),
       col=   c("black", "orange", "blue", "red"), lty=1:2, cex=0.8,
       title="Line types", text.font=4)
m
