function outputs = OCP_meltingspeed_Ipopt_pc(input,filename)

% Input data
nrun = input.nrun;
nc = input.nc;
endtime = input.endtime_Itf;
Tb_def = input.Tb;  % default heater temperature 
Tb_min = input.Tmin;  % lower bound
Tb_max = input.Tmax;  % upper bound
tb = linspace(0,input.endtime_Itf,nc+1)';
twall = zeros(nrun,1);
Tb = zeros(nc,nrun);
fval = zeros(nrun,1);

% Initial guess
if isempty(filename)
    Tb_ini = Tb_def*ones(nc,1);
    twall0 = 0;
else
    Data0 = load(filename).Data;
    nc0 = height(Data0.Tb);
    Tb0 = Data0.Tb;
    twall0 = Data0.twall;
    tb0 = linspace(0,endtime,nc0+1)';
    Tb_ini = interp1(tb0 , Tb0, tb);
end

% Optimization using Ipopt
for k = 1:nrun

    tic
    Tb_lb = Tb_min*ones(nc,1);
    Tb_ub = Tb_max*ones(nc,1);
    obj = @(Tb) obj_meltingspeed_pc(Tb, tb, input.dSdt_obj, input.IC, input.option_ode2, input);
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