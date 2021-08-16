
# this file contains the code for the tolerance study of lsoda for the 
# time problem with and without event detection

library(deSolve)
library(ggplot2)
library(plyr)

tolerances <- c(1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-11)

## defining the required functions
### no event detection
measures_implemented <- 0
direction <- "up"

model_no_event <- function(t, state, parms) {
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


experiment_no_event <- function(atol, rtol) {
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

    out <- ode(func=model_no_event, y=y0, times=tspan, parms=NULL, atol=atol, rtol=rtol)
    diag <- diagnostics(out)

    colnames(out) <- c("time",
                        paste("S", abs(round(log10(atol))), sep="_"),
                        paste("E", abs(round(log10(atol))), sep="_"),
                        paste("I", abs(round(log10(atol))), sep="_"),
                        paste("R", abs(round(log10(atol))), sep="_")

    )
    list(out=out, nfev=( diag$istate[3]), nstep=(diag$istate[2]))
}

### getting the data
res <- experiment_no_event(1e-7, 1e-7)
ans_no_event <- as.data.frame(res$out)
efficiency_no_event <- c(1e-7, res$nfev, res$nstep)
for (tolerance in tolerances) {
    res <- experiment_no_event(tolerance, tolerance)
    ans_no_event <- join(ans_no_event, as.data.frame(res$out),  by="time", type="left")
    row <- c(tolerance, res$nfev, res$nstep)
    efficiency_no_event <- rbind(efficiency_no_event, row)
}

### plotting the data
tolerances <- c(1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
colors <- c("1e-01" = "green", 
            "1e-02" = "blue", 
            "1e-04" = "orange",
            "1e-06" = "brown",
            "1e-07" = "yellow",
            "1e-08" = "red",
            "1e-10" = "purple",
            "1e-11" = "black")

plot_no_event <- ggplot(ans_no_event, aes(x=time)) + 
    geom_line(aes(y=E_1, color="1e-01"), size=1.2) + 
    geom_line(aes(y=E_2, color="1e-02"), size=1.2) + 
    geom_line(aes(y=E_4, color="1e-04"), size=1.2) +
    geom_line(aes(y=E_6, color="1e-06"), size=1.2)     +   
    geom_line(aes(y=E_7, color="1e-07"), size=1.2)     +   
    geom_line(aes(y=E_8, color="1e-08"), size=1.2) +  
    geom_line(aes(y=E_10, color="1e-10"), size=1.2) + 
    geom_line(aes(y=E_11, color="1e-11"), size=1.2) + 
    labs(x="time", y="E", color="legend") +
    scale_color_manual(values=colors)


#######################################################################

#### experiment with event  ####

library(deSolve)
library(ggplot2)
library(plyr)

tolerances <- c(1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-11)

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

experiment_with_event <- function(atol, rtol) {
    N <- 37.741 * (10^6)
    E0 <- 103
    I0 <- 1
    R0 <- 0
    S0 <- N - (E0 + I0 + R0)
    y0 <- c(S0, E0, I0, R0)
    t0 <- 0

    res <- NULL
    out <- NULL
    diag <- NULL
    sum_nfev <- 0
    sum_nstep <- 0

    res <- matrix(nrow=0, ncol=5)
    t_initial <- t0
    y_initial <- y0
    measures_implemented <- FALSE
    # need to stop at 179.5 so that radau can get a non-zero tspan
    while (t_initial < 179) {
        tspan = seq(from=t_initial, to=180, by=1)
        if (measures_implemented) {
            out <- ode(func=model_with_measures, y=y_initial, 
                    times=tspan, parms=NULL, rootfun=root_10000, 
                    atol=atol, rtol=rtol)
            diag <- diagnostics(out)
            measures_implemented <- FALSE
        } else {
            out <- ode(func=model_no_measures, y=y_initial, 
                   times=tspan, parms=NULL, rootfun=root_25000, 
                   atol=atol, rtol=rtol)
            diag <- diagnostics(out)
            measures_implemented <- TRUE
        }

        t_change <- length(out[, 1])
        t_initial <- (out[t_change, 1])
        y_initial <- out[t_change, ]
        y_initial <- y_initial[!(names(y_initial) %in% c("time"))]

        sum_nfev <- sum_nfev + diag$istate[3]
        sum_nstep <- sum_nstep + diag$istate[2]
        res <- rbind(res, out)
    }

    colnames(res) <- c("time", 
                    paste("S", abs(round(log10(atol))), sep="_"),
                    paste("E", abs(round(log10(atol))), sep="_"),
                    paste("I", abs(round(log10(atol))), sep="_"),
                    paste("R", abs(round(log10(atol))), sep="_")
                  )
    list(out=res, nfev=(sum_nfev), nstep=(sum_nstep))
}


tolerances_string <- c("1e-7", "1e-1", "1e-2", "1e-4", "1e-6", "1e-8", "1e-10", "1e-11")
colors <- c("green", "blue",  "orange", "brown", "yellow", "red", "purple", "black")
res <- experiment_with_event(1e-7, 1e-7)
length <- c(nrow(res$out))
i <- 1
plot(res$out[, 1], res$out[, 3], type="l", col=colors[i], lwd=2)
i <- i + 1
efficiency_event <- c(1e-7, res$nfev, res$nstep)
for (tolerance in tolerances) {
    res <- experiment_with_event(tolerance, tolerance)
    lines(res$out[, 1], res$out[, 3], col=colors[i], lwd=2)
    i <- i + 1
    row <- c(tolerance, res$nfev, res$nstep)
    efficiency_event <- rbind(efficiency_event, row)
}
legend(-6, 25000, 
   legend=tolerances_string,
   col=colors, lty=1:2, cex=0.8,
   title="Line types", text.font=4)

# printing efficiency data
efficiency_data <- cbind(efficiency_no_event, efficiency_event) # used to report efficiency
efficiency_data
