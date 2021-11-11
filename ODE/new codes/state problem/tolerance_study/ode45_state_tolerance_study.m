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
tolerances = [1e-1, 1e-2, 1e-4, 1e-6, 1e-7, 1e-9, 1e-10, 1e-11];

% ode45 without event detection
figure;
tolerance = tolerances(1);
measures_implemented = 0;
direction = up;
nfev = 0;
options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
[t,y] = ode45(@(t, y) f_with_if(t, y),tspan,y0, options);
plot(t, y(:, 2), 'lineWidth', 1.5);
hold on;
nfevs = [nfev];

for i = 2:length(tolerances)
    tolerance = tolerances(i);
    measures_implemented = 0;
    direction = up;
    nfev = 0;
    options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
    [t,y] = ode45(@(t, y) f_with_if(t, y), tspan, y0, options);
    plot(t, y(:, 2), 'lineWidth', 1.5);
    hold on;
    nfevs(i) = nfev;
end
hold off;
legend("1e-1", "1e-2", "1e-4", "1e-6", "1e-7", "1e-9", "1e-10", "1e-11");
A = [tolerances; nfevs];
transpose(A)

% with event detection
t0 = 0;
t_final = 180;

% ode45 experiment with event detection
figure;
tolerance = tolerances(1);
res = [S0 E0 I0 R0];
t_res = [t0];
y_initial = y0;
t_initial = t0;
nfev = 0;
measures_implemented = 0;
while t_initial < t_final
    tspan = [t_initial t_final];
    if (measures_implemented == 0)
        opts = odeset('Events', @g_25000, 'RelTol', tolerance, 'AbsTol', tolerance);
        [t,y,te,ye,ie] = ode45(@(t, y) f_no_measures(t, y), tspan, y_initial, opts);
        measures_implemented = 1;
    else 
        opts = odeset('Events', @g_10000, 'RelTol', tolerance, 'AbsTol', tolerance);
        [t,y,te,ye,ie] = ode45(@(t, y) f_measures(t, y), tspan, y_initial, opts);
        measures_implemented = 0;
    end 
    y_initial = y(end, :);
    t_initial = t(end);
    t_res = vertcat(t_res, t);
    res = vertcat(res, y);
end
plot(t_res, res(:, 2), 'lineWidth', 1.5);
hold on;
nfevs = [nfev];


for i = 2:length(tolerances)
    tolerance = tolerances(i);

    res = [S0 E0 I0 R0];
    t_res = [t0];
    y_initial = y0;
    t_initial = t0;
    nfev = 0;
    measures_implemented = 0;
    while t_initial < t_final
        tspan = [t_initial t_final];
        if (measures_implemented == 0)
            opts = odeset('Events', @g_25000, 'RelTol', tolerance, 'AbsTol', tolerance);
            [t,y,te,ye,ie] = ode45(@(t, y) f_no_measures(t, y), tspan, y_initial, opts);
            measures_implemented = 1;
        else 
            opts = odeset('Events', @g_10000, 'RelTol', tolerance, 'AbsTol', tolerance);
            [t,y,te,ye,ie] = ode45(@(t, y) f_measures(t, y), tspan, y_initial, opts);
            measures_implemented = 0;
        end 
        y_initial = y(end, :);
        t_initial = t(end);
        t_res = vertcat(t_res, t);
        res = vertcat(res, y);
    end
    plot(t_res, res(:, 2), 'lineWidth', 1.5);
    hold on;
    nfevs(i) = nfev;
end
hold off;
legend("1e-1", "1e-2", "1e-4", "1e-6", "1e-7", "1e-9", "1e-10", "1e-11");
A = [tolerances; nfevs];
transpose(A)

% functions 
function res = f_with_if(t, y)
    global measures_implemented direction nfev up down;
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
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

    alpha = 1.0/8.0;
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
    dRdt = gamma*I - mu*R;
    res = [dSdt; dEdt; dIdt; dRdt];
end

% event detection function
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
