global nfev;
nfev = 0;
N = 37.741d6;

tspan = [0 95];
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

nfev = 0;
[t_ode45,y_ode45] = ode45(@(t, y) f_with_if(t, y),tspan,y0);
count_ode45 = nfev;

nfev = 0;
[t_ode15s,y_ode15s] = ode15s(@(t, y) f_with_if(t, y),tspan,y0);
count_ode15s = nfev;

% setting the tolerance
% options = odeset('RelTol',1e-8,'AbsTol',1e-10);
% [t,y] = ode45(odefun,tspan,y0,options)

plot(t_ode45, y_ode45(:, 3), t_ode15s, y_ode15s(:, 3), 'lineWidth', 1.5);
legend("ode45", "ode15s");
count_ode45
count_ode15s

function res = f_with_if(t, y)
    global nfev ;
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    N = 37.741d6;
    t_event = 27;
    alpha = 1.0/8.0;
    beta = 0.9;
    if (t > t_event) 
        beta = 0.005;
    end
    gamma = 0.06;
    mu = 0.01/365;

    dSdt = mu*N - mu*S - (beta/N)*I*S;
    dEdt = (beta/N)*I*S - alpha*E - mu*E;
    dIdt = alpha*E - gamma*I - mu*I;
    dRdt = gamma*I - mu*R; 
    res = [dSdt; dEdt; dIdt; dRdt];
end
