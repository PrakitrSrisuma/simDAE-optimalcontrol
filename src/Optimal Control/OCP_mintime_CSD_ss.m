function outputs = OCP_mintime_CSD_ss(input)

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

% Formulate the NLP
Xk = input.IC;
for k=0:nc-1
    % New NLP variable for the control
    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, Tmin];
    ubw = [ubw,  Tmax];
    w0 = [w0,  Tb_ini];
    Fk = F('x0',Xk,'p', Uk);

    % Integrate till the end of the interval
    Xk = Fk.xf(:,end);
    S = Fk.xf(end,:);
    J=J+S(end)^2;

end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}));
solver = nlpsol('solver', 'ipopt', prob, input.optCSD);

% Solve the NLP
tic; sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw); twall = toc;
w_opt = full(sol.x);

% Export data
Data.Tb_opt = [w_opt;w_opt];
Data.tb = [0;tb];
Data.twall = twall;
Data.ext = '';
Data.Tb_ini = w0;
outputs = Data;

return