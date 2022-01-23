

// CAN ONLY DO UP TO 1e-7

tolerances = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];
colors = ["red", "blue", "green", "yellow", "pink", "purple", "black", "orange"];
tolerances_string = ["1e-1", "1e-2", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7"];
up = 1;
down = 0;

global("measures_implemented", "direction", "nfev")
measures_implemented = 0;
direction = up;
nfev = 0;
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
tspan = [0:1:90];

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
    nfev = 0
    y_rkf = ode("rkf", y0, t0, tspan, atol, rtol, model_with_if);
    E_rkf = y_rkf(2, :);
    plot(tspan, E_rkf, "color", color);
    A = [atol, nfev]
end
legend(tolerances_string);
xlabel("time", "fontsize", 4);
ylabel("E(t)", "fontsize", 4);