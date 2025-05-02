function outputs = OCP_tempchange_Ipopt(input,filename)

% Input data
nrun = input.nrun;
nc = input.nc;
endtime = input.endtime_temp;
tspan = (0:input.dt_plot:endtime)';
Tb_min = input.Tmin;
Tb_max = input.Tmax;
Tb_def = 0.5*(Tb_max+Tb_min);
tb = linspace(0,input.endtime_Itf,nc+1)';
twall = zeros(nrun,1);
Tb = zeros(nc+1,nrun);
fval = zeros(nrun,1);

% Initial guess
if isempty(filename)
    Tb_ini = Tb_def*ones(nc+1,1);
    twall0 = 0;
else
    Data0 = load(filename).Data;
    nc0 = height(Data0.Tb_opt)-1;
    Tb0 = Data0.Tb_opt;
    twall0 = Data0.twall;
    tb0 = linspace(0,endtime,nc0+1)';
    Tb_ini = interp1(tb0 , Tb0, tb);
end

% Optimization using Ipopt
for k = 1:nrun
    tic
    Tb_lb = Tb_min*ones(nc+1,1);
    Tb_ub = Tb_max*ones(nc+1,1);
    obj = @(Tb) obj_tempchange_pl(Tb, tb, tspan, input.dTdt_obj, input.IC, input.option_ode2, input);
    Optimizer_Ipopt = opti('obj',obj,'lb',Tb_lb,'ub',Tb_ub,'options',input.optIpopt);
    [Tb(:,k), fval(k)] = solve(Optimizer_Ipopt,Tb_ini);
    twall(k) = toc;
end

% Export data
Data.Tb_opt = Tb(:,1);
Data.tb = tb;
Data.twall = mean(twall) + twall0;
Data.ext = filename;
Data.Tb_ini = Tb_ini;
outputs = Data;

return