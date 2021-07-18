
global("nfev")
nfev=0

function res = f_no_measures(t, y)
    global("nfev")
    nfev = nfev + 1

    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    N = 37.741d6;  
    alpha = 1.0/8.0;
    beta = 0.9;
    gamma = 0.06;
    mu = 0.01/365;

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    res(3) = alpha*E - gamma*I - mu*I;
    res(4) = gamma*I - mu*R;
endfunction

function res = f_measures(t, y)
    global("nfev")
    nfev = nfev + 1
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    N = 37.741d6;      
    alpha = 1.0/8.0;
    beta = 0.005;
    gamma = 0.06;
    mu = 0.01/365;

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    res(3) = alpha*E - gamma*I - mu*I;
    res(4) = gamma*I - mu*R;
endfunction

function sol = g_25000(t, y)
    E = y(2);
    sol = E - 25000;
endfunction

function sol = g_10000(t, y)
    E = y(2);
    sol = E - 10000;
endfunction

N = 37.741d6; 
ng = 1;
t0 = 0;
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];
t_final = 180;

res = [1 1 1 1];
y_initial = y0;
t_initial = t0;
measures_implemented = 0;
while t_initial < t_final
    tspan = [t_initial:1:t_final];
    if (measures_implemented == 0)
        [y, rd] = ode("root", y_initial, t_initial, tspan, f_no_measures, ng, g_25000);
        measures_implemented = 1;
    else 
        [y, rd] = ode("root", y_initial, t_initial, tspan, f_measures, ng, g_10000);
        measures_implemented = 0;
    end 
    y = y';
    t_change = length(y(:, 1));
    t_initial = t_initial + t_change;
    y_initial = y(t_change, :)';
    res = cat(1, res, y);
end

xset("thickness",2);
E = res(:, 2);
t_span = (0:t_final);
plot(t_span, E);

nfev

// ode_root USES lsodar 
