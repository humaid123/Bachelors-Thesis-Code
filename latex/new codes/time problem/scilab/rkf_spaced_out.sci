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
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

tspan = [0:1:95];
nfev = 0;
y_lsoda = ode(y0, t0, tspan, f_with_if);
I_lsoda = y_lsoda(3, :);
count_lsoda = nfev;

tspan_step_1 = [0:1:95];
nfev = 0;
y_rkf_step_1 = ode("rkf", y0, t0, tspan_step_1, f_with_if);
I_rkf_step_1 = y_rkf_step_1(3, :);
count_rkf_step_1 = nfev;

tspan_step_2 = [0:2:95];
nfev = 0;
y_rkf_step_2 = ode("rkf", y0, t0, tspan_step_2, f_with_if);
I_rkf_step_2 = y_rkf_step_2(3, :);
count_rkf_step_2 = nfev;

tspan_step_5 = [0:5:95];
nfev = 0;
y_rkf_step_5 = ode("rkf", y0, t0, tspan_step_5, f_with_if);
I_rkf_step_5 = y_rkf_step_5(3, :);
count_rkf_step_5 = nfev;

tspan_step_7 = [0:7:95];
nfev = 0;
y_rkf_step_7 = ode("rkf", y0, t0, tspan_step_7, f_with_if);
I_rkf_step_7 = y_rkf_step_7(3, :);
count_rkf_step_7 = nfev;

scf(0);
xset("thickness", 2);
plot(tspan, I_lsoda, "color", "red");
plot(tspan_step_1, I_rkf_step_1, "color", "brown");
plot(tspan_step_2, I_rkf_step_2, "color", "green");
plot(tspan_step_5, I_rkf_step_5, "color", "blue");
plot(tspan_step_7, I_rkf_step_7, "color", "black");
legend(["ANS", "step=1", "step=2", "step=5", "step=7"]);
A = [count_rkf_step_1, count_rkf_step_2, count_rkf_step_5, count_rkf_step_7]

// ############################################################
clear;

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

E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0_before = [S0; E0; I0; R0];

t0_before = 0;
tspan_before = [0:1:27];
t0_after = 27;
tspan_after = [27:1:95];
tspan = cat(2, tspan_before, tspan_after);

// lsoda
nfev = 0;
y_lsoda_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_lsoda_before(:, length(tspan_before));
y_lsoda_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_lsoda = cat(2, y_lsoda_before, y_lsoda_after);
I_lsoda = y_lsoda(3, :);
count_lsoda = nfev;


// rkf tstep 1
t0_before = 0;
tspan_before = [0:1:27];
t0_after = 27;
tspan_after = [27:1:95];
tspan_step_1 = cat(2, tspan_before, tspan_after);

nfev = 0;
y_tstep1_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_tstep1_before(:, length(tspan_before));
y_tstep1_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_tstep1 = cat(2, y_tstep1_before, y_tstep1_after);
I_tstep1 = y_tstep1(3, :);
count_tstep1 = nfev;

// rkf tstep 2
t0_before = 0;
tspan_before = [0:2:27, 27];
t0_after = 27;
tspan_after = [27:2:95];
tspan_step_2 = cat(2, tspan_before, tspan_after);

nfev = 0;
y_tstep2_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_tstep2_before(:, length(tspan_before));
y_tstep2_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_tstep2 = cat(2, y_tstep2_before, y_tstep2_after);
I_tstep2 = y_tstep2(3, :);
count_tstep2 = nfev;

// rkf tstep 5
t0_before = 0;
tspan_before = [0:5:27, 27];
t0_after = 27;
tspan_after = [27:5:95];
tspan_step_5 = cat(2, tspan_before, tspan_after);

nfev = 0;
y_tstep5_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_tstep5_before(:, length(tspan_before));
y_tstep5_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_tstep5 = cat(2, y_tstep5_before, y_tstep5_after);
I_tstep5 = y_tstep5(3, :);
count_tstep5 = nfev;

// rkf tstep 7
t0_before = 0;
tspan_before = [0:7:27, 27];
t0_after = 27;
tspan_after = [27:7:95];
tspan_step_7 = cat(2, tspan_before, tspan_after);

nfev = 0;
y_tstep7_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_tstep7_before(:, length(tspan_before));
y_tstep7_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_tstep7 = cat(2, y_tstep7_before, y_tstep7_after);
I_tstep7 = y_tstep7(3, :);
count_tstep7 = nfev;

scf(1);
xset("thickness", 2);
plot(tspan, I_lsoda, "color", "red");
plot(tspan_step_1, I_tstep1, "color", "brown");
plot(tspan_step_2, I_tstep2, "color", "green");
plot(tspan_step_5, I_tstep5, "color", "blue");
plot(tspan_step_7, I_tstep7, "color", "black");
legend(["ANS", "step=1", "step=2", "step=5", "step=7"]);
A = [count_tstep1, count_tstep2, count_tstep5, count_tstep7]
