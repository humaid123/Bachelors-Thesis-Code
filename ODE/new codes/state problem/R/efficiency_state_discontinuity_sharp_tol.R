
library(deSolve)
library(ggplot2)
library(plyr)

atol <- 1e-13
rtol <- 1e-13

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
out <- ode(func=model_with_if, y=y0, times=tspan, parms=NULL, atol=atol, rtol=rtol)
toc <- Sys.time()
res <- diagnostics(out)
count <- c("lsoda", as.integer(res$istate[3]), toc-tic)
m <- count

colnames(out) <- c("time", "S_lsoda", "E_lsoda", "I_lsoda", "R_lsoda")
final <- as.data.frame(out)

methods <- c("daspk", "euler", "rk4", "ode45", "radau", "bdf", "adams")
for (method in methods) {
    .GlobalEnv$measures_implemented <- 0
    .GlobalEnv$direction <- "up"
    tic <- Sys.time()
    out <- ode(func=model_with_if, y=y0, times=tspan, parms=NULL, atol=atol, rtol=rtol, method=method)
    toc <- Sys.time()
    res <- diagnostics(out)
    count <- c(method, as.integer(res$istate[3]), toc-tic)
    m <- rbind(m, count)

    colnames(out) <- c("time",
                        paste("S", method, sep="_"),
                        paste("E", method, sep="_"),
                        paste("I", method, sep="_"),
                        paste("R", method, sep="_")
    )

    final <- join(as.data.frame(out), final, by="time", type="left") 
}

colors <- c("lsoda" = "green", 
            "daspk" = "yellow", 
            "euler" = "orange",
            "rk4"   = "blue",
            "ode45" = "pink",
            "radau" = "yellow",
            "bdf"   = "brown",
            "adams" = "red")

p <- ggplot(final, aes(x=time)) + 
    geom_line(aes(y=E_lsoda, color="lsoda"), size=1.2) +  
    geom_line(aes(y=E_daspk, color="daspk"), size=1.2) +  
    geom_line(aes(y=E_euler, color="euler"), size=1.2) +
    geom_line(aes(y=E_rk4, color="rk4"), size=1.2)     +   
    geom_line(aes(y=E_ode45, color="ode45"), size=1.2) +  
    geom_line(aes(y=E_bdf, color="bdf"), size=1.2)     +  
    geom_line(aes(y=E_adams, color="adams"), size=1.2) +  
    labs(x="time", y="E", color="legend") +
    scale_color_manual(values=colors) 
p
m

h <- ggplot(final, aes(x=time)) + 
    geom_line(aes(y=E_radau, color="radau"), size=1.2) 
h
