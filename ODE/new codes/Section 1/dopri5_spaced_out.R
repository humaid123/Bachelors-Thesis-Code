
library(deSolve)
nfev <- 0
model <- function(t, state, parms) {
    S <- state[1]
    E <- state[2]
    I <- state[3]
    R <- state[4]

    nfev <<- nfev + 1
    
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
nfev <- 0
lsoda <- ode(func=model, y=y0, times=tspan, parms=NULL)
nfev


# ode45 spaced out
nfev <- 0
tspan <- seq(from=0, to=95, by=1)
tstep_1 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_1)
count <- c("tstep_1", as.integer(res$istate[3]), nfev)
m <- count

nfev <- 0
tspan <- seq(from=0, to=95, by=3)
tstep_3 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_3)
count <- c("tstep_3", as.integer(res$istate[3]), nfev)
m <- rbind(m, count)

nfev <- 0
tspan <- seq(from=0, to=95, by=5)
tstep_5 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_5)
count <- c("tstep_5", as.integer(res$istate[3]), nfev)
m <- rbind(m, count)

nfev <- 0
tspan <- seq(from=0, to=95, by=7)
tstep_7 <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, parms=NULL)
res <- diagnostics(tstep_7)
count <- c("tstep_7", as.integer(res$istate[3]), nfev)
m <- rbind(m, count)


plot(tstep_7[, 1], tstep_7[, 4], type="l", col="black", lwd=2, 
    xlab="time", 
    ylab="I(t)", 
    cex.lab = 1.5,
    cex.axis = 1.5,
    cex.main = 1.5,
    cex.sub = 1.5
)
lines(tstep_5[, 1], tstep_5[, 4], type="l", col="brown", lwd=2)
lines(tstep_3[, 1], tstep_3[, 4], col="orange", lwd=2)
lines(tstep_1[, 1], tstep_1[, 4], col="blue", lwd=2)
lines(lsoda[, 1], lsoda[, 4], type="l", col="red", lwd=2)
legend(1, 30000, 
       legend=c("space=7", "space=5", "space=3", "space=1", "ANS"),
       col=   c("black", "brown", "orange", "blue", "red"), lty=1, lwd=2,
       title="Legend", text.font=4, cex=1.1)
m

# Notes on plot() in R
# the xlab and ylab in plot() gives the labels of the axes
# lty=1 says that the  legend contains solid lines, lty=2 => dashed lines
# lwd says the width of the lines in the legends



# # ode45 spaced out with an hmax
# nfev <- 0
# tspan <- seq(from=0, to=95, by=1)
# tstep_1_hmax <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, hmax=1e60, parms=NULL)
# res <- diagnostics(tstep_1_hmax)
# count <- c("tstep_1_with_hmax", as.integer(res$istate[3]), nfev)
# m <- rbind(m, count)

# nfev <- 0
# tspan <- seq(from=0, to=95, by=3)
# tstep_3_hamx <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, hmax=1e60, parms=NULL)
# res <- diagnostics(tstep_3_hmax)
# count <- c("tstep_3_with_hmax", as.integer(res$istate[3]), nfev)
# m <- rbind(m, count)

# nfev <- 0
# tspan <- seq(from=0, to=95, by=5)
# tstep_5_hmax <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, hmax=1e60, parms=NULL)
# res <- diagnostics(tstep_5_hmax)
# count <- c("tstep_5_with_hmax", as.integer(res$istate[3]), nfev)
# m <- rbind(m, count)

# nfev <- 0
# tspan <- seq(from=0, to=95, by=7)
# tstep_7_hmax <- ode(func=model, y=y0, times=tspan, method="ode45", atol=1e-1, rtol=1e-1, hmax=1e60, parms=NULL)
# res <- diagnostics(tstep_7_hmax)
# count <- c("tstep_7_with_hmax", as.integer(res$istate[3]), nfev)
# m <- rbind(m, count)



# plot(tstep_7_hmax[, 1], tstep_7_hmax[, 4], type="l", col="black", lwd=2)
# lines(lsoda[, 1], lsoda[, 4], type="l", col="red", lwd=2)
# lines(tstep_5_hmax[, 1], tstep_5_hmax[, 4], type="l", col="brown", lwd=2)
# lines(tstep_3_hmax[, 1], tstep_3_hmax[, 4], col="orange", lwd=2)
# lines(tstep_1_hmax[, 1], tstep_1_hmax[, 4], col="blue", lwd=2)
# legend(1, 30000, 
#        legend=c("hmax_space=7", "hmax_space=5", "hmax_space=3", "hmax_space=1", "ANS"),
#        col=   c("black", "orange", "blue", "red"), lty=1:2, cex=0.8,
#        title="Line types", text.font=4)
# m