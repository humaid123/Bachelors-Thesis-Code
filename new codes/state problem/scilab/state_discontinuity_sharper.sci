
up = 1;
down = 0;

global("measures_implemented", "direction", "nfev")
measures_implemented = 0;
direction = up;
nfev=0

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

t0 = 0;
tspan = [0:2:180];

E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];

atol = 1e-13
rtol = 1e-13

// lsoda
nfev = 0;
y_lsoda = ode(y0, t0, tspan, atol, rtol, model_with_if);
E_lsoda = y_lsoda(2, :);
count_lsoda = nfev

// stiff
nfev = 0
y_stiff = ode("stiff", y0, t0, tspan, atol, rtol, model_with_if);
E_stiff = y_stiff(2, :);
count_stiff = nfev

// rkf - DOES NOT ALLOW SHARP TOLERANCES
// y_rkf = ode("rkf", y0, t0, tspan, atol, rtol, model_with_if);
// E_rkf = y_rkf(2, :);
// E_rkf(91) = E_rkf(90); // required as rkf messes up one last iteration

// rk
nfev = 0
y_rk = ode("rk", y0, t0, tspan, atol, rtol, model_with_if);
E_rk = y_rk(2, :);
count_rk = nfev

// adams
nfev = 0
y_adams = ode("adams", y0, t0, tspan, atol, rtol, model_with_if);
E_adams = y_adams(2, :);
count_adams = nfev

xset("thickness",2);
A = [E_lsoda; E_stiff; E_rk; E_adams]';
plot2d(tspan, A, [1, 2, 3, 4, 5]);
h1=legend(["lsoda", "stiff", "rk", "adams"]);
[count_lsoda; count_stiff; count_rk; count_adams]
// rkf CANNOT BE USED
