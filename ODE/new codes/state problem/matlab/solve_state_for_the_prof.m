global nfev;
nfev = 0;

N = 37.741d6; 
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];
t0 = 0;
t_final = 95;

% ode45 experiment
res = [S0 E0 I0 R0];
t_res = [t0];
y_initial = y0;
t_initial = t0;
nfev = 0;
measures_implemented = 0;
while t_initial < t_final
    tspan = [t_initial t_final];
    if (measures_implemented == 0)
        opts = odeset('Events', @g_25000);
        [t,y,te,ye,ie] = ode45(@f_no_measures,tspan,y_initial,opts);
        measures_implemented = 1;
    else 
        opts = odeset('Events', @g_10000);
        [t,y,te,ye,ie] = ode45(@f_measures,tspan,y_initial,opts);
        measures_implemented = 0;
    end 
    y_initial = y(end, :);
    t_initial = t(end);
    t_res = vertcat(t_res, t);
    res = vertcat(res, y);
end
t_ode45 = t_res;
y_ode45 = res;
nfev_ode45 = nfev;

% ode15s experiment
res = [S0 E0 I0 R0];
t_res = [t0];
y_initial = y0;
t_initial = t0;
nfev = 0;
measures_implemented = 0;
while t_initial < t_final
    tspan = [t_initial t_final];
    if (measures_implemented == 0)
        opts = odeset('Events', @g_25000);
        [t,y,te,ye,ie] = ode15s(@f_no_measures,tspan,y_initial,opts);
        measures_implemented = 1;
    else 
        opts = odeset('Events', @g_10000);
        [t,y,te,ye,ie] = ode15s(@f_measures,tspan,y_initial,opts);
        measures_implemented = 0;
    end 
    y_initial = y(end, :);
    t_initial = t(end);

    t_res = vertcat(t_res, t);
    res = vertcat(res, y);
end
t_ode15s = t_res;
y_ode15s = res;
nfev_ode15s = nfev;

% plot(t_ode45, y_ode45(:, 2), t_ode15s, y_ode15s(:, 2), 'lineWidth', 1.5);
% only ode45
plot(t_ode45, y_ode45(:, 2), 'lineWidth', 1.5);
%legend("ode45", "ode15s");
xlabel("time");
ylabel("E(t)");
nfev_ode45
nfev_ode15s

function res = f_no_measures(t, y)
    global nfev;
    nfev = nfev + 1;

    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    N = 37.741d6;  
    alpha = 1.0/8.0;
    beta = 0.9;
    gamma = 0.06;
    mu = 0.01/365;

    dSdt = mu*N - mu*S - (beta/N)*I*S;
    dEdt = (beta/N)*I*S - alpha*E - mu*E;
    dIdt = alpha*E - gamma*I - mu*I;
    dRdt = gamma*I - mu*R;
    res = [dSdt; dEdt; dIdt; dRdt];
end

function res = f_measures(t, y)
    global nfev;
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    N = 37.741d6;      
    alpha = 1.0/8.0;
    beta = 0.005;
    gamma = 0.06;
    mu = 0.01/365;

    dSdt = mu*N - mu*S - (beta/N)*I*S;
    dEdt = (beta/N)*I*S - alpha*E - mu*E;
    dIdt = alpha*E - gamma*I - mu*I;
    dRdt = gamma*I - mu*R;
    res = [dSdt; dEdt; dIdt; dRdt];
end

function [value,isterminal,direction] = g_25000(t, y)
    E = y(2);
    value = E - 25000;
    isterminal = 1;
    direction = 0; 
end

function [value,isterminal,direction] = g_10000(t, y)
    E = y(2);
    value = E - 10000;
    isterminal = 1;
    direction = 0; 
end
