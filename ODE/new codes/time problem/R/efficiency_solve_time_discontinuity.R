
library(deSolve)
library(ggplot2)
library(plyr)

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

t0 <- 0
tspan_before <- seq(from=0, to=27, by=1)
tspan_after <- seq(from=27, to=95, by=1)

N <- 37.741 * (10^6)
E0 <- 103
I0 <- 1
R0 <- 0
S0 <- N - (E0 + I0 + R0)
y0_before <- c(S0, E0, I0, R0)

# lsoda is the default function
# need to do this outside the loop to get the initial final
# data frame
t1 <- Sys.time()
out_before <- ode(func=model_before, y=y0_before, times=tspan_before, parms=NULL)
t2 <- Sys.time()
y0_after <- (out_before[nrow(out_before), ])[!(colnames(out_before) %in% c("time"))]
t3 <- Sys.time()
out_after <- ode(func=model_after, y=y0_after, times=tspan_after, parms=NULL)
t4 <- Sys.time()

res_before <- diagnostics(out_before)
res_after <- diagnostics(out_after)

count <- c("lsoda", res_before$istate[3] + res_after$istate[3], t2-t1 + t4-t3)
m <- count
out <- rbind(out_before, out_after)
colnames(out) <- c("time", "S_lsoda", "E_lsoda", "I_lsoda", "R_lsoda")
final <- as.data.frame(out)



##### OTHER METHODS ####################

experiment <- function(method, dense) {
    t0 <- 0
    tspan_before <- seq(from=0, to=27, by=1)
    tspan_after <- seq(from=27, to=95, by=1)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)

    y0_before <- c(S0, E0, I0, R0)

    t1 <- Sys.time()
    out_before <- ode(func=model_before, y=y0_before, 
                          times=tspan_before, parms=NULL, method=method)
    t2 <- Sys.time()
    y0_after <- (out_before[nrow(out_before),])[!(colnames(out_before) %in% c("time"))]
    t3 <- Sys.time()

    out_after <- ode(func=model_after, y=y0_after, 
                         times=tspan_after, parms=NULL, method=method)
    t4 <- Sys.time()

    res_before <- diagnostics(out_before)
    res_after <- diagnostics(out_after)
    count <- c(method, res_before$istate[3] + res_after$istate[3], t2-t1 + t4-t3)

    out <- rbind(out_before, out_after)
    colnames(out) <- c("time",
                        paste("S", method, sep="_"),
                        paste("E", method, sep="_"),
                        paste("I", method, sep="_"),
                        paste("R", method, sep="_")
    )
    list(out=out, count=count, timeElapsed=(t2-t1 + t4-t3))
}

methods <- c("daspk", "euler", "rk4", "ode45", "radau", "bdf", "adams")
for (method in methods) {
    res <- experiment(method, FALSE)
    out <- res$out
    count <- res$count
    timeElapsed <- res$timeElapsed
    m <- rbind(m, count)
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
    geom_line(aes(y=I_lsoda, color="lsoda"), size=1.2) +  # lsoda
    geom_line(aes(y=I_daspk, color="daspk"), size=1.2) +  # on top of lsoda
    geom_line(aes(y=I_euler, color="euler"), size=1.2) +
    geom_line(aes(y=I_rk4, color="rk4"), size=1.2)     +  # on top of lsoda 
    geom_line(aes(y=I_ode45, color="ode45"), size=1.2) +  # on top of lsoda
    geom_line(aes(y=I_radau, color="radau"), size=1.2) +  # on top of lsoda
    geom_line(aes(y=I_bdf, color="bdf"), size=1.2)     +  # on top of lsoda
    geom_line(aes(y=I_adams, color="adams"), size=1.2) +  # on top of lsoda
    labs(x="time", y="I(t)", color="legend") +
    scale_color_manual(values=colors)  + 
    theme(text = element_text(size = 20)) 
p
m
