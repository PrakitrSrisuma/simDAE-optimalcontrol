function outputs = OCP_mintime_Ipopt(input)

% Input data
nrun = input.nrun;
endtime = 100;
tspan = [0 endtime];
Tb_min = input.Tmin; 
Tb_max = input.Tmax;  
Tb_def = 0.5*(Tb_max+Tb_min);  
twall = zeros(nrun,1);
Tb = zeros(1,nrun);
fval = zeros(nrun,1);
options = input.option_ode2;
Tb_ini = Tb_def;

% Optimization using Ipopt
for k = 1:nrun
    tic
    Tb_lb = Tb_min;
    Tb_ub = Tb_max;
    obj = @(Tb) obj_mintime(Tb, tspan, input.IC, options, input);
    Optimizer_Ipopt = opti('obj',obj,'lb',Tb_lb,'ub',Tb_ub,'options',input.optIpopt);
    [Tb(:,k), fval(k)] = solve(Optimizer_Ipopt,Tb_ini);
    twall(k) = toc;
end

% Export data
Data.Tb_opt = [Tb(1);Tb(1)];
[t,T] = ode15s(@(t,T)PDE_2Phases(t,T,Tb(1),input), tspan, input.IC, options);
Data.tb = [0;t(end)] ;
Data.twall = mean(twall);
Data.ext = '';
Data.Tb_ini = Tb_ini;
outputs = Data;

return