
global("nfev");
nfev = 0;

tolerances = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
colors = ["red", "blue", "green", "yellow", "black", "purple", "pink"];
tolerances_string = ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"];

N = 37.741d6;
t_event = 27;

function res = model_with_if(t, y)
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

function res = model_before(t, y)
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
    res(4) = gamma*I - mu*R;
endfunction

function res = model_after(t, y)
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

t0 = 0;
tspan = [0:1:95];
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

scf(0);
xset("thickness", 2);
for i = 1:length(tolerances)
    atol = tolerances(i);
    rtol = tolerances(i);
    color = colors(i);

    nfev = 0;
    y_rkf = ode("rkf", y0, t0, tspan, atol, rtol, model_with_if);
    I_rkf = y_rkf(3, :);
    plot(tspan, I_rkf, "color", color);
    A = [atol, nfev]
end
legend(tolerances_string);

scf(1);
t0_before = 0;
tspan_before = [0:1:27];
t0_after = 27;
tspan_after = [27:1:95];
tspan = cat(2, tspan_before, tspan_after);

y0_before = [S0; E0; I0; R0];
xset("thickness", 2);

for i = 1:length(tolerances)
    atol = tolerances(i);
    rtol = tolerances(i);
    color = colors(i);
    
    nfev = 0;
    y_rkf_before = ode("rkf", y0_before, t0_before, tspan_before, atol, rtol, model_before);
    y0_after = y_rkf_before(:, length(tspan_before));
    y_rkf_after = ode("rkf", y0_after, t0_after, tspan_after, atol, rtol, model_after);
    
    y_rkf = cat(2, y_rkf_before, y_rkf_after);
    I_rkf = y_rkf(3, :);
    plot(tspan, I_rkf, "color", color);
    A = [atol, nfev]
end
legend(tolerances_string);
