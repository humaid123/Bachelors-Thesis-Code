global nfev;
nfev = 0;
N = 37.741d6;

tspan_before = [0 27];
tspan_after = [27 95];

E0 = 103d0;
I0 = 1d0;
R0 = 0d0;
S0 = N - (E0 + I0 + R0);
y0_before = [S0; E0; I0; R0];

nfev = 0;
[t_ode45_before,y_ode45_before] = ode45(@(t, y) f_before(t, y),tspan_before,y0_before);
y0_after = y_ode45_before(end, :);
[t_ode45_after,y_ode45_after] = ode45(@(t, y) f_after(t, y),tspan_after,y0_after);
t_ode45 = vertcat(t_ode45_before, t_ode45_after);
y_ode45 = vertcat(y_ode45_before, y_ode45_after);
count_ode45 = nfev;

nfev = 0;
[t_ode15s_before,y_ode15s_before] = ode15s(@(t, y) f_before(t, y),tspan_before,y0_before);
y0_after = y_ode15s_before(end, :);
[t_ode15s_after,y_ode15s_after] = ode15s(@(t, y) f_after(t, y),tspan_after,y0_after);
t_ode15s = vertcat(t_ode15s_before, t_ode15s_after);
y_ode15s = vertcat(y_ode15s_before, y_ode15s_after);
count_ode15s = nfev;

plot(t_ode45, y_ode45(:, 3), t_ode15s, y_ode15s(:, 3), 'lineWidth', 1.5);
legend("ode45", "ode15s");
count_ode45
count_ode15s

function res = f_before(t, y)
    global nfev ;
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    beta = 0.9;
    gamma = 0.06;
    mu = 0.01/365;
    N = 37.741d6;

    dSdt = mu*N - mu*S - (beta/N)*I*S;
    dEdt = (beta/N)*I*S - alpha*E - mu*E;
    dIdt = alpha*E - gamma*I - mu*I;
    dRdt = gamma*I - mu*R; 
    res = [dSdt; dEdt; dIdt; dRdt];
end

function res = f_after(t, y)
    global nfev ;
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    beta = 0.005;
    gamma = 0.06;
    mu = 0.01/365;
    N = 37.741d6;

    dSdt = mu*N - mu*S - (beta/N)*I*S;
    dEdt = (beta/N)*I*S - alpha*E - mu*E;
    dIdt = alpha*E - gamma*I - mu*I;
    dRdt = gamma*I - mu*R; 
    res = [dSdt; dEdt; dIdt; dRdt];
end
