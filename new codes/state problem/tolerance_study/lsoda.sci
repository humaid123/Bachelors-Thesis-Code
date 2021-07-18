tolerances = [1e-1, 1e-2, 1e-4, 1e-6, 1e-7, 1e-8, 1e-10, 1e-11];
colors = ["red", "blue", "green", "yellow", "pink", "purple", "black", "orange"];
tolerances_string = ["1e-1", "1e-2", "1e-4", "1e-6", "1e-7", "1e-9", "1e-10", "1e-11"];


up = 1;
down = 0;

global("measures_implemented", "direction", "nfev")
measures_implemented = 0;
direction = up;
nfev=0;
N = 37.741d6;

function res = model_with_if(t, y)
    global('measures_implemented', 'direction', 'nfev');
    nfev = nfev + 1
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);
    
    alpha = 1.0/8.0;
    //disp(isglobal('direction'), direction, direction == up)
    if (direction == up) then
        if (E > 25000) then
            measures_implemented = 1;
            direction = down;
            //disp("switch to down", E, measures_implemented, direction)
        end
    elseif (direction == down) then
        if (E < 10000) then
            measures_implemented = 0;
            direction = up;
            //disp("switch to up", E, measures_implemented, direction)
        end
    end

    beta = 0.9
    if (measures_implemented == 1) 
        beta = 0.005
    end
    gamma = 0.06;
    mu = 0.01/365;
    //disp(t, y, beta, measures_implemented, direction)

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    res(3) = alpha*E - gamma*I - mu*I;
    res(4) = gamma*I - mu*R ;
endfunction

function res = f_no_measures(t, y)
    global('measures_implemented', 'direction', 'nfev');
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
    global('measures_implemented', 'direction', 'nfev');
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


t0 = 0;
tspan = [0:2:180];

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
    y_lsoda = ode(y0, t0, tspan, atol, rtol, model_with_if);
    E_lsoda = y_lsoda(2, :);
    plot(tspan, E_lsoda, "color", color);
    A = [atol, nfev]
end
legend(tolerances_string);

// USING LSODAR
ng = 1;
t_final = 180;

scf(1);
xset("thickness", 2);
for i = 1:length(tolerances)
    atol = tolerances(i);
    rtol = tolerances(i);
    color = colors(i);

    nfev = 0
    res = [1 1 1 1];
    y_initial = y0;
    t_initial = t0;
    measures_implemented = 0;
    while t_initial < t_final
        tspan = [t_initial:1:t_final];
        if (measures_implemented == 0)
            [y, rd] = ode("root", y_initial, t_initial, tspan, atol, rtol, f_no_measures, ng, g_25000);
            measures_implemented = 1;
        else 
            [y, rd] = ode("root", y_initial, t_initial, tspan, atol, rtol, f_measures, ng, g_10000);
            measures_implemented = 0;
        end 
        y = y';
        t_change = length(y(:, 1));
        t_initial = t_initial + t_change;
        y_initial = y(t_change, :)';
        res = cat(1, res, y);
    end

    E_lsoda = res(:, 2);
    t_span = (0:length(E_lsoda)-1);
    plot(t_span, E_lsoda, "color", color);
    A = [atol, nfev]
end
legend(tolerances_string);
