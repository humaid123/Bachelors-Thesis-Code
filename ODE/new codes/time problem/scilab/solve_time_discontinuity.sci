
global("nfev");
nfev = 0;
N = 37.741d6;

function res = f_before_event(t, y)
    global("nfev");
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    beta = 0.9;
    gamma = 0.06;
    mu = 0.01/365;

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    res(3) = alpha*E - gamma*I - mu*I;
    res(4) = gamma*I - mu*R ;
endfunction

function res = f_after_event(t, y)
    global("nfev");
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    beta = 0.005;
    gamma = 0.06;
    mu = 0.01/365;

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    res(3) = alpha*E - gamma*I - mu*I;
    res(4) = gamma*I - mu*R;
endfunction

t0_before = 0;
tspan_before = [0:1:27];
t0_after = 27;
tspan_after = [27:1:95];
tspan = cat(2, tspan_before, tspan_after);

E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0_before = [S0; E0; I0; R0];

// lsoda
nfev = 0;
y_lsoda_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_lsoda_before(:, length(tspan_before));
y_lsoda_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_lsoda = cat(2, y_lsoda_before, y_lsoda_after);
I_lsoda = y_lsoda(3, :);
count_lsoda = nfev;

// stiff
nfev = 0;
y_stiff_before = ode("stiff", y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_stiff_before(:, length(tspan_before));
y_stiff_after = ode("stiff", y0_after, t0_after, tspan_after, f_after_event);
y_stiff = cat(2, y_stiff_before, y_stiff_after);
I_stiff = y_stiff(3, :);
count_stiff = nfev;

// rkf
nfev = 0;
y_rkf_before = ode("rkf", y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_rkf_before(:, length(tspan_before));
y_rkf_after = ode("rkf", y0_after, t0_after, tspan_after, f_after_event);
y_rkf = cat(2, y_rkf_before, y_rkf_after);
I_rkf = y_rkf(3, :);
count_rkf = nfev;

// rk
nfev = 0;
y_rk_before = ode("rk", y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_rk_before(:, length(tspan_before));
y_rk_after = ode("rk", y0_after, t0_after, tspan_after, f_after_event);
y_rk = cat(2, y_rk_before, y_rk_after);
I_rk = y_rk(3, :);
count_rk = nfev;

// adams
nfev = 0;
y_adams_before = ode("adams", y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_adams_before(:, length(tspan_before));
y_adams_after = ode("adams", y0_after, t0_after, tspan_after, f_after_event);
y_adams = cat(2, y_adams_before, y_adams_after);
I_adams = y_adams(3, :);
count_adams = nfev;

xset("thickness",2);
A = [I_lsoda; I_stiff; I_rkf; I_rk; I_adams]';
plot2d(tspan, A, [1, 2, 3, 4, 5])
h1=legend(["lsoda", "stiff", "rkf", "rk", "adams"])
A = [count_lsoda; count_stiff; count_rkf; count_rk; count_adams]
