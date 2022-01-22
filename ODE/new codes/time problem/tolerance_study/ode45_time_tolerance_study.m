global nfev;
nfev = 0;
N = 37.741d6;
format shortG;

tspan = [0 95];
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

% setting the tolerance
tolerances = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
figure;

tolerance = tolerances(1);
nfev = 0;
options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
[t,y] = ode45(@(t, y) f_with_if(t, y), tspan, y0, options);
plot(t, y(:, 3), 'lineWidth', 1.5);
hold on;
nfevs = [nfev];

for i = 2:length(tolerances)
    tolerance = tolerances(i);
    nfev = 0;
    options = odeset('RelTol', tolerance, 'AbsTol', tolerance);
    [t,y] = ode45(@(t, y) f_with_if(t, y), tspan, y0, options);
    plot(t, y(:, 3), 'lineWidth', 1.5);
    hold on;
    nfevs(i) = nfev;
end
hold off;
legend("1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7");
xlabel("time");
ylabel("I(t)");
A = [tolerances; nfevs];
transpose(A)

% event tolerance study
tspan_before = [0 27];
tspan_after = [27 95];

E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0_before = [S0; E0; I0; R0];

% setting the tolerance
tolerances = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
figure;

tolerance = tolerances(1);
nfev = 0;
options = odeset('RelTol',tolerance,'AbsTol',tolerance);
[t_ode45_before,y_ode45_before] = ode45(@(t, y) f_before(t, y),tspan_before,y0_before,options);
y0_after = y_ode45_before(end, :);
[t_ode45_after,y_ode45_after] = ode45(@(t, y) f_after(t, y),tspan_after,y0_after,options);
t = vertcat(t_ode45_before, t_ode45_after);
y = vertcat(y_ode45_before, y_ode45_after);
plot(t, y(:, 3), 'lineWidth', 1.5);
hold on;
nfevs = [nfev];

for i = 2:length(tolerances)
    tolerance = tolerances(i);
    nfev = 0;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [t_ode45_before,y_ode45_before] = ode45(@(t, y) f_before(t, y),tspan_before,y0_before,options);
    y0_after = y_ode45_before(end, :);
    [t_ode45_after,y_ode45_after] = ode45(@(t, y) f_after(t, y),tspan_after,y0_after,options);
    t = vertcat(t_ode45_before, t_ode45_after);
    y = vertcat(y_ode45_before, y_ode45_after);
    plot(t, y(:, 3), 'lineWidth', 1.5);
    hold on;
    nfevs(i) = nfev;
end
hold off;
legend("1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7");
xlabel("time");
ylabel("I(t)");

A = [tolerances; nfevs];
transpose(A)

% function
function res = f_with_if(t, y)
    global nfev ;
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    t_event = 27;
    alpha = 1.0/8.0;
    beta = 0.9;
    if (t > t_event) 
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
