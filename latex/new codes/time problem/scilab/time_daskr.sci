
global("nfev");
nfev = 0;
N = 37.741d6;
t_event = 27;

function [res, ires] = f_with_if(t, y, ydot)
    global("nfev");
    nfev = nfev + 1;
    S = y(1);
    E = y(2);
    I = y(3);
    R = y(4);

    t_event = 27;    
    N = 37.741d6;
    alpha = 1.0/8.0;
    beta = 0.9;
    if (t > t_event) 
        beta = 0.005;
    end
    gamma = 0.06;
    mu = 0.01/365;

    disp([S, E, I, R])

    disp([t_event, N, alpha, beta, t, gamma, mu])

    res(1) = mu*N - mu*S - (beta/N)*I*S;
    disp(1, res(1))
    res(2) = (beta/N)*I*S - alpha*E - mu*E;
    disp(2, res(2))
    res(3) = alpha*E - gamma*I - mu*I;
    disp(3, res(3))
    res(4) = gamma*I - mu*R; 
    disp(4, res(4))

    disp(res)

    // if we reach here the estimation was done correctly
    ires = 0

    disp(ires)
endfunction

function c = surface(t, y)
    c = 1
endfunction

t0 = 0;
tspan = [0:1:95];
E0 = 103;
I0 = 1;
R0 = 0;
S0 = N - (E0 + I0 + R0);
y0 = [S0; E0; I0; R0];
// taken from the lsoda solution at time=1
y1 = [37740889.6138738; 97.0055598659178; 12.9546334006217; 0.425932906336397];
ydot0 = y1 - y0;
ng = 1;

info = list([], 0, [], [], [], 0, [], 0, [], 0, 0, [], [], 1);


// info(7)
// if ydot0 is set so that g(t0, y0, ydot0) = 0 then set 
// info(7)=[]. Otherwise, set info(7)=[+-1, ..., +-1], with 
// info(7)(i) = 1 if y(i) is a differential variable and 
// info(7)(i) = -1 if y(i) is an algebraic variable 
// (if its derivatives do not appear explicitly in the 
// function g(t, y, ydot)).
info=list([],0,[],[],[],0,0);
info(7)=1;


// dassl in red
//nfev = 0;
//[y_dae, nn] = dassl([y0, ydot0], t0, tspan, f_with_if, info) // ng, surf, info);
//I_dae = y_dae(3, :);
//count_dae = nfev;

//y = dae(y0, t0, tspan, f_with_if)
[y, rd] = dae("root", y0, t0, tspan, f_with_if, ng, surface)


//function [r, ires]=chemres(t, y, yd)
//    r(1) = -0.04*y(1) + 1d4*y(2)*y(3) - yd(1);
//    r(2) =  0.04*y(1) - 1d4*y(2)*y(3) - 3d7*y(2)*y(2) - yd(2);
//    r(3) =       y(1) +     y(2)      + y(3)-1;
//    ires =  0;
//endfunction

//function pd=chemjac(x, y, yd, cj)
//    pd = [-0.04-cj , 1d4*y(3)               , 1d4*y(2);
//           0.04    ,-1d4*y(3)-2*3d7*y(2)-cj ,-1d4*y(2);
//           1       , 1                      , 1       ]
//endfunction

//x0 = [1; 0; 0];
//xd0 = [-0.04; 0.04; 0];
//t = [1.d-5:0.02:.4, 0.41:.1:4, 40, 400, 4000, 40000, 4d5, 4d6, 4d7, 4d8, 4d9, 4d10];

//y = dae(x0, 0, t, chemres); // Returns requested observation time points

