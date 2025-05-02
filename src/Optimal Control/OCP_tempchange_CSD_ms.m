function outputs = OCP_tempchange_CSD_ms(input,filename)

import casadi.*
n1 = input.n1;
dTdt_obj = input.dTdt_obj;
dt_plot = input.dt_plot;
T = input.endtime_temp;
Tmin = input.Tmin;
Tmax = input.Tmax;
Tb_def = 0.5*(Tmax+Tmin); 
nc = input.nc;
tspan = [(0:dt_plot:T/nc)';T/nc]; 
tspan = unique(tspan);
tb = linspace(0,T,nc+1)';

% Initial guess
if isempty(filename)
    Tb_ini = Tb_def*ones(nc+1,1);
    twall0 = 0;
else
    Data0 = load(filename).Data;
    nc0 = height(Data0.Tb_opt)-1;
    Tb0 = Data0.Tb_opt;
    twall0 = Data0.twall;
    tb0 = linspace(0,T,nc0+1)';
    Tb_ini = interp1(tb0 , Tb0, tb);
end

% Declare model variables
x = MX.sym('x',2*n1+1);
u = MX.sym('u',2);
t = MX.sym('t');
x0 = input.IC;

% Formulate discrete time dynamics
dae = struct('x',x,'p',u,'t',t,'ode',PDE_2Phases_CasADi(x,t,T/nc,u,input));
F = integrator('F', 'cvodes', dae, 0, tspan);

% Evaluate at a test point
Fk = F('x0',x0,'p',Tb_def*ones(2,1));
S_test = full(evalf(Fk.xf(end,:)));
plot(tspan,S_test)

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
T_ini = cell(nc,1);
T0 = x0;
for i = 1:nc
    [t_tmp,T_tmp] = ode15s(@(t,T)PDE_2Phases(t,T,[[Tb_ini(i);Tb_ini(i+1)],[tspan(1);tspan(end)]],input), tspan, T0, input.option_ode2);
    T_ini{i} = T_tmp(end,:);
    T0 = T_ini{i};
end

% Formulate the NLP
for k=0:nc-1

    % New NLP variable for the control
    if k == 0
        Uk = MX.sym(['U_' num2str(k)],2);
        w = {w{:}, Uk};
        lbw = [lbw; Tmin; Tmin];
        ubw = [ubw;  Tmax; Tmax];
        w0 = [w0;  Tb_ini(1); Tb_ini(2)];
        Fk = F('x0',Xk,'p', Uk);
        Uf = Uk(2);
    else
        Uk = MX.sym(['U_' num2str(k)]);
        w = {w{:}, Uk};
        lbw = [lbw; Tmin];
        ubw = [ubw;  Tmax];
        w0 = [w0;  Tb_ini(k+2)];
        Fk = F('x0',Xk,'p', [Uf;Uk]);
        Uf = Uk;
    end

    % Integrate till the end of the interval
    Xk_end = Fk.xf(:,end);
    sol = mean(Fk.xf(n1+1:end-1,:));
    J=J+obj_tempchange_CasADi(tspan,sol,dTdt_obj);

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
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
% solver = nlpsol('solver', 'ipopt', prob); 
% options.ipopt.tol = 1e-6;
solver = nlpsol('solver', 'ipopt', prob, input.optCSD); 

% Solve the NLP
tic; sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg); twall = toc;
w_opt = full(sol.x);

% Plot the solution
tb = linspace(0,input.endtime_Itf,nc+1)';
Tb_opt = zeros(nc+1,1);
for i = 1:nc+1
    if i == 1
        Tb_opt(i) = w_opt(i*(2*n1+1)+i);
    else
        Tb_opt(i) = w_opt((i-1)*(2*n1+1)+i);
    end
end

% Export data
Data.Tb_opt = Tb_opt;
Data.tb = tb;
Data.twall = twall + twall0;
Data.ext = filename;
Data.Tb_ini = w0;
outputs = Data;

return