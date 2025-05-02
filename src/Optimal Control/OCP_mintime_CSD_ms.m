function outputs = OCP_mintime_CSD_ms(input)

% Input data
import casadi.*
n1 = input.n1;
dt_plot = input.dt_plot;
T = 6.9;
Tmin = input.Tmin;
Tmax = input.Tmax;
Tb_def = 0.5*(Tmax+Tmin); 
nc = input.nc;
tspan = [(0:dt_plot:T/nc)';T/nc]; 
tspan = unique(tspan);
tb = linspace(0,T,nc)';
Tb_ini = Tb_def*ones(nc,1);

% Declare model variables
x = MX.sym('x',2*n1+1,1);
u = MX.sym('u',1);
t = MX.sym('t');
x0 = input.IC;

% CVODES from the SUNDIALS suite
dae = struct('x',x,'p',u,'ode',PDE_2Phases_CasADi_pc(x,u,input));
F = integrator('F', 'cvodes', dae, 0, tspan);

% Evaluate at a test point
Fk = F('x0',x0,'p',Tb_def);
S_test = full(evalf(Fk.xf(end,:)));

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
Xk = MX.sym('X0',2*n1+1);
w = {w{:}, Xk};
lbw = [lbw; x0];
ubw = [ubw; x0];
w0 = [w0; x0];

% Initial guess for states
T0 = x0;
[t_tmp,T_tmp] = ode15s(@(t,T)PDE_2Phases(t,T,Tb_ini,input), tspan, T0, input.option_ode2);
T_ini{1} = T_tmp(end,:);

% Formulate the NLP
Xk = input.IC;
for k=0:nc-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw; Tmin];
    ubw = [ubw;  Tmax];
    w0 = [w0;  Tb_ini];
    Fk = F('x0',Xk,'p', Uk);

    % Integrate till the end of the interval
    Xk_end = Fk.xf(:,end);
    S = Fk.xf(end,:);
    J=J+S(end)^2;

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 2*n1+1);
    w = [w, {Xk}];
    lbw = [lbw; -inf*ones(2*n1,1); input.tol];
    ubw = [ubw;  inf*ones(2*n1,1); 1];
    w0 = [w0; T_ini{k+1}'];

    % Add equality constraint
    g = [g, {Xk_end-Xk}];
    lbg = [lbg; zeros(2*n1+1,1)];
    ubg = [ubg; zeros(2*n1+1,1)];

end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}));
solver = nlpsol('solver', 'ipopt', prob, input.optCSD);

% Solve the NLP
tic; sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw); twall = toc;
w_opt = full(sol.x);
w_opt = w_opt(2*n1+2);

% Export data
Data.Tb_opt = [w_opt;w_opt];
Data.tb = [0;tb];
Data.twall = twall;
Data.ext = '';
Data.Tb_ini = w0;
outputs = Data;

return