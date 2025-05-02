function outputs = OCP_tempchange_CSD_ss(input,filename)

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
x = MX.sym('x',2*n1+1,1);
u = MX.sym('u',2);
t = MX.sym('t');
x0 = input.IC;

% CVODES from the SUNDIALS suite
dae = struct('x',x,'p',u,'t',t,'ode',PDE_2Phases_CasADi(x,t,T/nc,u,input));
F = integrator('F', 'cvodes', dae, 0, tspan);

% Evaluate at a test point
Fk = F('x0',x0,'p',Tb_def*ones(2,1));
Temp = full(evalf(Fk.xf(n1+1:end-1,:)));

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
    if k == 0
        Uk = MX.sym(['U_' num2str(k)],2);
        w = {w{:}, Uk};
        lbw = [lbw, Tmin, Tmin];
        ubw = [ubw,  Tmax, Tmax];
        w0 = [w0,  Tb_ini(1), Tb_ini(2)];
        Fk = F('x0',Xk,'p', Uk);
        Uf = Uk(2);
    else
        Uk = MX.sym(['U_' num2str(k)]);
        w = {w{:}, Uk};
        lbw = [lbw, Tmin];
        ubw = [ubw,  Tmax];
        w0 = [w0,  Tb_ini(k+2)];
        Fk = F('x0',Xk,'p', [Uf;Uk]);
        Uf = Uk;
    end

    % Integrate till the end of the interval
    Xk = Fk.xf(:,end);
    sol = mean(Fk.xf(n1+1:end-1,:));
    J=J+obj_tempchange_CasADi(tspan,sol,dTdt_obj);

end

% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}));
solver = nlpsol('solver', 'ipopt', prob, input.optCSD);

% Solve the NLP
tic; sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw); twall = toc;
w_opt = full(sol.x);

% Export data
Data.Tb_opt = w_opt;
Data.tb = tb;
Data.twall = twall + twall0;
Data.ext = filename;
Data.Tb_ini = w0;
outputs = Data;

return