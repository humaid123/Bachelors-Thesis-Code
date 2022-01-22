
# this file contains the code for the tolerance study of lsoda for the 
# time problem with and without event detection

library(deSolve)
library(ggplot2)
library(plyr)

tolerances <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)


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


experiment_no_event <- function(atol, rtol) {
    t0 <- 0
    tspan <- seq(from=0, to=95, by=1)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0 <- c(S0, E0, I0, R0)

    sum_nfev <- 0
    sum_nstep <- 0
    sum_time <- 0
    num_iter <- 10
    out <- NULL
    for (i in as.list(seq(from=1, to=num_iter, by=1))) {
        t1 <- Sys.time()
        out <- ode(func=model_no_event, y=y0, times=tspan, parms=NULL, 
                   atol=atol, rtol=rtol)
        t2 <- Sys.time()
        diag <- diagnostics(out)
        sum_nfev <- sum_nfev + diag$istate[3]
        sum_nstep <- sum_nstep + diag$istate[2]
        sum_time <- sum_time + t2-t1
    }
    colnames(out) <- c("time",
                        paste("S", abs(round(log10(atol))), sep="_"),
                        paste("E", abs(round(log10(atol))), sep="_"),
                        paste("I", abs(round(log10(atol))), sep="_"),
                        paste("R", abs(round(log10(atol))), sep="_")

    )
    list(out=out, nfev=(sum_nfev/num_iter), nstep=(sum_nstep/num_iter), timeElapsed=(sum_time/num_iter))
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

experiment_with_event <- function(atol, rtol) {
    t0 <- 0
    tspan_before <- seq(from=0, to=27, by=1)
    tspan_after <- seq(from=27, to=95, by=1)
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0_before <- c(S0, E0, I0, R0)

    sum_time <- 0
    sum_nfev <- 0
    sum_nstep <- 0
    num_iter <- 10
    out_before <- NULL
    out_after <- NULL
    for (i in as.list(seq(from=1, to=num_iter, by=1))) {
        t1 <- Sys.time()
        out_before <- ode(func=model_before, y=y0_before, times=tspan_before, 
                          parms=NULL, atol=atol, rtol=rtol)
        t2 <- Sys.time()

        y0_after <- (out_before[nrow(out_before),])[!(colnames(out_before) %in% c("time"))]
        t3 <- Sys.time()
        out_after <- ode(func=model_after, y=y0_after, times=tspan_after, 
                         parms=NULL, atol=atol, rtol=rtol)
        t4 <- Sys.time()

        diag_before <- diagnostics(out_before)
        diag_after <- diagnostics(out_after)
        sum_nfev <- sum_nfev + (diag_before$istate[3] + diag_after$istate[3])
        sum_nstep <- sum_nstep + (diag_before$istate[2] + diag_after$istate[2])
        sum_time <- sum_time + (t2-t1 + t4-t3)
    }
    out <- rbind(out_before, out_after)
    colnames(out) <- c("time",
                        paste("S", abs(round(log10(atol))), sep="_"),
                        paste("E", abs(round(log10(atol))), sep="_"),
                        paste("I", abs(round(log10(atol))), sep="_"),
                        paste("R", abs(round(log10(atol))), sep="_")

    )
    list(out=out, nfev=(sum_nfev/num_iter), nstep=(sum_nstep/num_iter), timeElapsed=(sum_time/num_iter))
}

### getting the data
res <- experiment_no_event(1e-7, 1e-7)
ans_no_event <- as.data.frame(res$out)
efficiency_no_event <- c(1e-7, res$nfev, res$timeElapsed, res$nstep)
for (tolerance in tolerances) {
    res <- experiment_no_event(tolerance, tolerance)
    ans_no_event <- join(as.data.frame(res$out), ans_no_event, by="time", type="left")
    row <- c(tolerance, res$nfev, res$timeElapsed, res$nstep)
    efficiency_no_event <- rbind(efficiency_no_event, row)
}

res <- experiment_with_event(1e-7, 1e-7)
ans_event <- as.data.frame(res$out)
efficiency_event <- c(1e-7, res$nfev, res$timeElapsed, res$nstep)
for (tolerance in tolerances) {
    res <- experiment_with_event(tolerance, tolerance)
    ans_event <- join(as.data.frame(res$out), ans_event, by="time", type="left")
    row <- c(tolerance, res$nfev, res$timeElapsed, res$nstep)
    efficiency_event <- rbind(efficiency_event, row)
}

### plotting the data
colors <- c("1e-01" = "green", 
            "1e-02" = "blue", 
            "1e-03" = "orange",
            "1e-04"   = "yellow",
            "1e-05" = "red",
            "1e-06" = "black",
            "1e-07"   = "brown")

plot_no_event <- ggplot(ans_no_event, aes(x=time)) + 
    geom_line(aes(y=I_1, color="1e-01"), size=1.2) + 
    geom_line(aes(y=I_2, color="1e-02"), size=1.2) + 
    geom_line(aes(y=I_3, color="1e-03"), size=1.2) +
    geom_line(aes(y=I_4, color="1e-04"), size=1.2)     +   
    geom_line(aes(y=I_5, color="1e-05"), size=1.2) +  
    geom_line(aes(y=I_6, color="1e-06"), size=1.2) + 
    geom_line(aes(y=I_7, color="1e-07"), size=1.2)  + 
    labs(x="time", y="I(t)", color="legend") +
    scale_color_manual(values=colors)  + 
    theme(text = element_text(size = 20))  

plot_with_event <- ggplot(ans_event, aes(x=time)) + 
    geom_line(aes(y=I_1, color="1e-01"), size=1.2) + 
    geom_line(aes(y=I_2, color="1e-02"), size=1.2) + 
    geom_line(aes(y=I_3, color="1e-03"), size=1.2) +
    geom_line(aes(y=I_4, color="1e-04"), size=1.2)     +   
    geom_line(aes(y=I_5, color="1e-05"), size=1.2) +  
    geom_line(aes(y=I_6, color="1e-06"), size=1.2) + 
    geom_line(aes(y=I_7, color="1e-07"), size=1.2)  + 
    labs(x="time", y="I(t)", color="legend") +
    scale_color_manual(values=colors)  + 
    theme(text = element_text(size = 20))  

# printing efficiency data
efficiency_data <- cbind(efficiency_no_event, efficiency_event) # used to report efficiency
efficiency_data
