
global("nfev");
nfev = 0;
N = 37.741d6;
t_event = 27;

function res = f_with_if(t, y)
    global("nfev");
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    beta = 0.9;
    if (t > t_event) 
        beta = 0.005;
    end
    gamma = 0.06;
    mu = 0.01/365;

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    res(3) = alpha*E - gamma*I - mu*I;
    res(4) = gamma*I - mu*R; 
endfunction

t0 = 0;
tspan = [0:1:95];
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

// lsoda in red
nfev = 0;
y_lsoda = ode(y0, t0, tspan, f_with_if);
I_lsoda = y_lsoda(3, :);
count_lsoda = nfev;

// stiff=bdf in green
nfev = 0;
y_stiff = ode("stiff", y0, t0, tspan, f_with_if);
I_stiff = y_stiff(3, :);
count_stiff = nfev;

// rkf = ode(45) in blue
nfev = 0;
y_rkf = ode("rkf", y0, t0, tspan, f_with_if);
I_rkf = y_rkf(3, :);
count_rkf = nfev;

// rk = classical order 4 rk => black
nfev = 0;
y_rk = ode("rk", y0, t0, tspan, f_with_if);
I_rk = y_rk(3, :);
count_rk = nfev;

// adams method in pink
nfev = 0;
y_adams = ode("adams", y0, t0, tspan, f_with_if);
I_adams = y_adams(3, :);
count_adams = nfev;

xset("thickness",2);
A = [I_lsoda; I_stiff; I_rkf; I_rk; I_adams]';
plot2d(tspan, A, [1, 2, 3, 4, 5]);
legends(["lsoda", "stiff", "rkf", "rk", "adams"], [1, 2, 3, 4, 5],font_size=4);
xlabel("time", "fontsize", 4)
ylabel("I(t)", "fontsize", 4)

A = [count_lsoda; count_stiff; count_rkf; count_rk; count_adams]
