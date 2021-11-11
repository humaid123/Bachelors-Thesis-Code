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

t0 = 0;
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

tspan_step_7 = [0:2:95];
nfev = 0;
y_rkf_step_7 = ode("rkf", y0, t0, tspan_step_7, f_with_if);
I_rkf_step_7 = y_rkf_step_7(3, :);
count_rkf_no_event = nfev;


t0_before = 0;
tspan_before = [0:2:27, 27];
t0_after = 27;
tspan_after = [27:2:95];
tspan_step_7 = cat(2, tspan_before, tspan_after);
nfev = 0;
y_tstep7_before = ode(y0_before, t0_before, tspan_before, f_before_event);
y0_after = y_tstep7_before(:, length(tspan_before));
y_tstep7_after = ode(y0_after, t0_after, tspan_after, f_after_event);
y_tstep7 = cat(2, y_tstep7_before, y_tstep7_after);
I_tstep7 = y_tstep7(3, :);
count_rkf_with_event = nfev;

A = [count_rkf_no_event, count_rkf_with_event]
