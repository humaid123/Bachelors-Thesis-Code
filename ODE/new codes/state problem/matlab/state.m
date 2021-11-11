global measures_implemented direction nfev up down;
up = 1;
down = 0;
measures_implemented = 0;
direction = up;
nfev = 0;

N = 37.741d6;
tspan = [0 180];
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

measures_implemented = 0;
direction = up;
nfev = 0;
[t_ode45, y_ode45] = ode45(@(t, y) f_with_if(t, y), tspan, y0);
count_ode45 = nfev;

measures_implemented = 0;
direction = up;
nfev = 0;
[t_ode15s,y_ode15s] = ode15s(@(t, y) f_with_if(t, y), tspan, y0);
count_ode15s = nfev;

plot(t_ode45, y_ode45(:, 2), t_ode15s, y_ode15s(:, 2), 'lineWidth', 1.5);
legend("ode45", "ode15s");
count_ode45
count_ode15s


function res = f_with_if(t, y)
    global measures_implemented direction nfev up down;
    nfev = nfev + 1;

    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    if (direction == up)
        if (E > 25000)
            measures_implemented = 1;
            direction = down;
        end
    elseif (direction == down)
        if (E < 10000)
            measures_implemented = 0;
            direction = up;
        end
    end

    beta = 0.9;
    if (measures_implemented == 1) 
        beta = 0.005;
    end
    gamma = 0.06;
    mu = 0.01/365;
    N = 37.741d6;

    dSdt = mu*N - mu*S - (beta/N)*I*S;
    dEdt = (beta/N)*I*S - alpha*E - mu*E;
    dIdt = alpha*E - gamma*I - mu*I;
    dRdt = gamma*I - mu*R ;
    res = [dSdt; dEdt; dIdt; dRdt];
end
