% ==============================================================================
% Simulation-based Approach for Fast Optimal Control of a Stefan Problem 
% with Application to Cell Therapy
%
% This is a top-level routine.
% This routine considers an optimal control problem that controls the interface
% velocity in cell thawing, corresponding to Problem 3 in the paper.
% ==============================================================================
close all; clear; clc;

%% Preparation
% Add paths
addpath('Python','Events','Input Data','Plotting','Objective Functions', ...
    'Optimal Control','PDEs','Calculations')
addpath(genpath('CasADi'))  % add your path to CasADi here
addpath(genpath('OPTI-master'))  % add your path to IPOPT here

% Typical input data
ip0 = get_inputdata;  % unpack all default input data
ip0.nc = 12;  % control interval
ip = input_processing(ip0);  % processing the input  
n1 = ip.n1;  % number of grid points
Tb_def = ip.Tb;  % default heater temperature 
endtime = ip.endtime_Itf;  % final time
tspan = unique([(0:ip.dt_plot:endtime)'; endtime]);  % time span with interval = dt_plot
T_ini = ip.IC;  % initial condition
dSdt_obj = ip.dSdt_obj;  % target dSdt
S_target = 1+tspan*dSdt_obj;  % target interface position

% ODE solvers
tol_ode = ip.tol_ode;
option_ode1 = ip.option_ode1;
option_ode2 = ip.option_ode2;

% Optimal control
nrun = ip.nrun;  % number of simulation runs
nc = ip.nc;  % number of control intervals

%% Option selection ('on' or 'off')
orgsim = 'off';  % original simulation
Ipopt = 'off';  % Ipopt
fmincon = 'off';  % fmincon with a piecewise linear control
fmincon_par = 'off';  % fmincon with a piecewise linear control via parallel computing
CasADi_ss = 'off';  % CasADi single shooting
CasADi_ms = 'off';  % CasADi multiple shooting
DAE = 'off';  % simulation-based approach using GEKKO's solver called from Python
DAE2 = 'off';  % simulation-based approach with active upper bound
GEKKO = 'off';  % optimization with GEKKO
Ipopt_pc = 'off';  % Ipopt with a piecewise constant control
all = 'on';  % run all simulations (use to test all functions)
none = 'off';  % run nothing (use to debug code in the prep block)

switch all
case 'on'
    [orgsim,Ipopt,fmincon,fmincon_par,CasADi_ss,CasADi_ms,DAE,DAE2,GEKKO,Ipopt_pc] = deal('on');
end

switch none
case 'on'
    [orgsim,Ipopt,fmincon,fmincon_par,CasADi_ss,CasADi_ms,DAE,DAE2,GEKKO,Ipopt_pc] = deal('off');
end

%% Original simulation (no control)
switch orgsim
case 'on'
[t_sim,T_sim] = ode15s(@(t,T)PDE_2Phases(t,T,0,ip), tspan, T_ini, option_ode2);
S_sim = T_sim(:,end);  % interface position
dTdt = cal_dTdt(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)));
dSdt = cal_dSdt(t_sim,S_sim);

figure; plot_interface(t_sim,S_sim)
figure; plot_Tavg(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)))

end

%% Optimization using IPOPT with a piecewise linear control
switch Ipopt
case 'on'
disp('Solving the optimal control problem with Ipopt')

% Solve the optimal control problem
Data = OCP_meltingspeed_Ipopt(ip,'');

% Solve the ODEs with the optimal profiles
[t_Ipopt,T_Ipopt] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_Ipopt = T_Ipopt(:,end);
Data.RMSE = sqrt(sum((S_Ipopt-S_target).^2)/length(t_Ipopt));
Data.e2 = cal_dSdt_2norm(t_Ipopt,S_Ipopt,dSdt_obj);
disp(['e2_Ipopt = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_Ipopt,S_Ipopt) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_Ipopt,S_Ipopt,ip)

end

%% Optimization using fmincon with a piecewise linear control
switch fmincon
case 'on'
disp('Solving the optimal control problem with fmincon')

% Solve the optimal control problem
Data = OCP_meltingspeed_fmincon(ip,'');

% Solve the ODEs with the optimal profiles
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Data.RMSE = sqrt(sum((S_fmincon-S_target).^2)/length(t_fmincon));
Data.e2 = cal_dSdt_2norm(t_fmincon,S_fmincon,dSdt_obj);
disp(['e2_fmincon = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_fmincon,S_fmincon,ip)

end

%% Optimization using fmincon with a piecewise linear control via parallel computing
switch fmincon_par
case 'on'
disp('Solving the optimal control problem with fmincon via parallel computing')

% Solve the optimal control problem
Data = OCP_meltingspeed_fmincon_par(ip,'');

% Solve the ODEs with the optimal profiles 
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Data.RMSE = sqrt(sum((S_fmincon-S_target).^2)/length(t_fmincon));
Data.e2 = cal_dSdt_2norm(t_fmincon,S_fmincon,dSdt_obj);
disp(['e2_fmincon = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_fmincon,S_fmincon,ip)

end

%% Direct single shooting with CasADi
switch CasADi_ss
case 'on'
disp('Solving the optimal control problem with CasADi, single shooting')

% Solve the optimal control problem
Data = OCP_meltingspeed_CSD_ss(ip,'');

% Solve the ODEs with the optimal profiles
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
S_target = 1+t_CSD*dSdt_obj;
Data.RMSE = sqrt(sum((S_CSD-S_target).^2)/length(t_CSD));
Data.e2 = cal_dSdt_2norm(t_CSD,S_CSD,ip.dSdt_obj);

end

%% Direct multiple shooting with CasADi
switch CasADi_ms
case 'on'
disp('Solving the optimal control problem with CasADi, multiple shooting')

% Solve the optimal control problem
Data = OCP_meltingspeed_CSD_ms(ip,'');

% Solve the ODEs with the optimal profiles
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
S_target = 1+t_CSD*dSdt_obj;
Data.RMSE = sqrt(sum((S_CSD-S_target).^2)/length(t_CSD));
Data.e2 = cal_dSdt_2norm(t_CSD,S_CSD,ip.dSdt_obj);

end

%% Simulation-based approach
switch DAE
case 'on'
disp('Solving the optimal control problem with the simulation-based technique')

% Load the Python file
py_file = 'DAE_MeltingSpeed.py';

% DAE Solver from Python
Tb0 = 0.2186;  % consistent initial condition
cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
output_py = pyrunfile(py_file,'output_MATLAB',z1=ip.dSdt_obj,z2=Tb0,z3=ip.endtime_Itf,z4=3);
Theta_opt_DAE = double(output_py{2})';
Data.Tb_opt = ip.temp_dim(Theta_opt_DAE);
Data.tb = double(output_py{3})';
Data.twall = output_py{1};
Data.ext = '';
Data.Tb_ini = Tb0;
cd(fileparts(matlab.desktop.editor.getActiveFilename))
    
% Solve the ODEs with the optimal profiles
[t_DAE,T_DAE] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_DAE = T_DAE(:,end);
Data.RMSE = sqrt(sum((S_DAE-S_target).^2)/length(t_DAE));
Data.e2 = cal_dSdt_2norm(t_DAE,S_DAE ,dSdt_obj);
disp(['e2_DAE = ' , num2str(Data.e2)])

figure; plot_interface(t_DAE,S_DAE) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_DAE,S_DAE,ip)

end

%% Simulation-based approach with active upper bound
switch DAE2
case 'on'
disp('Solving the optimal control problem with the simulation-based technique')

% Load the Python file
py_file = 'DAE_MeltingSpeed_mixed.py';

% Pre-processing
ip2 = get_inputdata;
ip2.dSdt_obj = -0.15;
ip2.endtime_Itf = 7;
ip2 = input_processing(ip2);  % processing the input  
endtime = ip2.endtime_Itf;  % final time
tspan = [(0:ip2.dt_plot:endtime)'; endtime];  % time span with interval = dt_plot
tspan = unique(tspan);  % remove duplicated value
tspan_mixed = linspace(0,endtime,20);
dt_mixed = tspan_mixed(2)-tspan_mixed(1);
T_ini = ip2.IC;  % initial condition
dSdt_obj2 = ip2.dSdt_obj;  % target dSdt
S_target2 = 1+tspan*dSdt_obj2;  % target interface position
tb_master = [];
Tb_master = [];
twall_master = [];
Tb0 = 0.327581631457616;  % inconsistent initial condition
tcheck = 0;
count = 0;
nodes = 3;

% DAE Solver from Python
tic
while tcheck < endtime
    cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
    output_py = pyrunfile(py_file,'output_MATLAB',z1=dSdt_obj2,z2=Tb0,z3=tspan_mixed,z4=nodes,z5=T_ini);
    Theta_opt_DAE = double(output_py{2})';
    Tb_opt = ip2.temp_dim(Theta_opt_DAE);
    tb = double(output_py{3})';
    twall = output_py{1};

    sw = 0;
    
    for j = 2:length(Theta_opt_DAE)-1
        if Theta_opt_DAE(j)>=1 && Theta_opt_DAE(j+1)<1
            sw = j;
        end
    end

    if sw == 0
        sw = length(tb);
    end

    if count > 0
        Tb_master = [Tb_master;Theta_opt_DAE(2:sw)];
        Tb_master(Tb_master>1) = 1;
        tb_master = [tb_master;tb(2:sw)];
    else
        Tb_master = [Tb_master;Theta_opt_DAE(1:sw)];
        Tb_master(Tb_master>1) = 1;
        tb_master = [tb_master;tb(1:sw)];
    end
    cd(fileparts(matlab.desktop.editor.getActiveFilename))

    [t_tmp,T_tmp] = ode15s(@(t,T)PDE_2Phases(t,T,[ip2.temp_dim(Tb_master),tb_master],ip2), tb_master, ip2.IC, option_ode2);
    T_ini = T_tmp(end,:);
    Tb0 = 1;
    tspan_mixed = [tb(sw);endtime];
    nodes = 3;
    count = count+1;
    tcheck = tb_master(end);
end
    
% Solve the ODEs with the optimal profiles
T_ini = ip2.IC;
[t_DAE,T_DAE] = ode15s(@(t,T)PDE_2Phases(t,T,[ip2.temp_dim(Tb_master),tb_master],ip2), tspan, T_ini, option_ode2);
S_DAE = T_DAE(:,end);
Data.Tb_opt = ip2.temp_dim(Tb_master);
Data.tb = tb_master;
Data.twall = toc;
Data.ext = '';
Data.Tb_ini = Tb0;
Data.RMSE = sqrt(sum((S_DAE-S_target2).^2)/length(t_DAE));
Data.e2 = cal_dSdt_2norm(t_DAE,S_DAE ,dSdt_obj2);
disp(['e2_DAE = ' , num2str(Data.e2)])

figure; plot_interface(t_DAE,S_DAE) 
figure; plot_Tb(tb_master,Tb_master)
figure; plot_dSdt(t_DAE,S_DAE,ip2)

end

%% Optimization using GEKKO with a piecewise linear control
switch GEKKO
case 'on'
disp('Solving the optimal control problem with GEKKO')

% Load the Python file
py_file = 'Gekko_MeltingSpeed.py';
Tb0 = [];
cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
output_py = pyrunfile(py_file,'output_MATLAB',x=dSdt_obj);
Theta_opt_Gekko = double(output_py{2})';
Data.Tb_opt = ip.temp_dim(Theta_opt_Gekko);
Data.tb = double(output_py{3})';
Data.twall = output_py{1};
Data.ext = '';
Data.Tb_ini = Tb0;
cd(fileparts(matlab.desktop.editor.getActiveFilename))
    
% Solve the ODEs with the optimal profiles
[t_Gekko,T_Gekko] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_Gekko = T_Gekko(:,end);
Data.RMSE = sqrt(sum((S_Gekko-S_target).^2)/length(t_Gekko));
Data.e2 = cal_dSdt_2norm(t_Gekko,S_Gekko,dSdt_obj);
disp(['e2_Gekko = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_Gekko,S_Gekko) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dSdt(t_Gekko,S_Gekko,ip)

end

%% Optimization using IPOPT with a piecewise constant control
switch Ipopt_pc
case 'on'

% Solve the optimal control problem
Data = OCP_meltingspeed_Ipopt_pc(ip,'');
tb = Data.tb;
Tb = Data.Tb_opt;
t_Ipopt = [];
S_Ipopt = [];

for i = 1:ip.nc
    tspan = (tb(i):ip.dt_plot:tb(i+1));
    [t,T] = ode15s(@(t,T)PDE_2Phases(t,T,Tb(i),ip), tspan, T_ini, option_ode2);
    if i~=ip.nc
        t_Ipopt = [t_Ipopt;t(1:end-1)];
        S_Ipopt = [S_Ipopt;T(1:end-1,end)];
    else
        t_Ipopt = [t_Ipopt;t];
        S_Ipopt = [S_Ipopt;T(:,end)];
    end
    T_ini = T(end,:);
end
S_target = 1+t_Ipopt*dSdt_obj;
Data.RMSE = sqrt(sum((S_Ipopt-S_target).^2)/length(t_Ipopt));
Data.e2 = cal_dSdt_2norm(t_Ipopt,S_Ipopt,ip.dSdt_obj);
disp(['e2_Ipopt = ' , num2str(Data.e2)])


% Plot the result
figure; plot_interface(t_Ipopt,S_Ipopt) 
figure; plot_Tb_stairs(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_abserror(t_Ipopt,S_Ipopt,S_target)

end