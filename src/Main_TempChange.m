% ==============================================================================
% Simulation-based Approach for Fast Optimal Control of a Stefan Problem 
% with Application to Cell Therapy
%
% This is a top-level routine.
% This routine considers an optimal control problem that controls the rate of
% temperature change in cell thawing, corresponding to Problem 2 in the paper.
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
ip0.nc = 16;  % control interval
ip = input_processing(ip0);  % processing the input  
Tb_def = ip.Tb;  % default heater temperature 
n1 = ip.n1;  % number of grid points
endtime = ip.endtime_temp;  % final time
tspan = unique([(0:ip.dt_plot:endtime)'; endtime]);  % time span with interval = dt_plot
T_ini = ip.IC;  % initial condition
dTdt_obj = ip.dTdt_obj;  % target dTdt
T_target = tspan*dTdt_obj;  % target average temperature

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
GEKKO = 'off';  % optimization with GEKKO
all = 'on';  % run all simulations (use to test all functions)
none = 'off';  % run nothing (use to debug code in the prep block)

switch all
case 'on'
    [orgsim,Ipopt,fmincon,fmincon_par,CasADi_ss,CasADi_ms,DAE,GEKKO] = deal('on');
end

switch none
case 'on'
    [orgsim,Ipopt,fmincon,fmincon_par,CasADi_ss,CasADi_ms,DAE,GEKKO] = deal('off');
end

%% Original simulation (no control)
switch orgsim
case 'on'
[t_sim,T_sim] = ode15s(@(t,T)PDE_2Phases(t,T,0,ip), tspan, T_ini, option_ode2);
S_sim = T_sim(:,end);  % interface position
dTdt = cal_dTdt(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)));

figure; plot_interface(t_sim,S_sim)
figure; plot_Tavg(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)))

end

%% Optimization using IPOPT with a piecewise linear control
switch Ipopt
case 'on'
disp('Solving the optimal control problem with Ipopt')

% Solve the optimal control problem
Data = OCP_tempchange_Ipopt(ip,'');

% Solve the ODEs with the optimal profiles 
[t_Ipopt,T_Ipopt] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_Ipopt = T_Ipopt(:,end);
Tavg = cal_Tavg(T_Ipopt(:,n1+1:end-1));
dTdt = cal_dTdt(t_Ipopt,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_Ipopt));
Data.e2 = cal_dTdt_2norm(t_Ipopt,Tavg ,dTdt_obj);
disp(['e2_Ipopt = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_Ipopt,S_Ipopt) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_Ipopt,cal_Tavg(T_Ipopt(:,n1+1:end-1)),ip)

end

%% Optimization using fmincon with a piecewise linear control
switch fmincon
case 'on'
disp('Solving the optimal control problem with fmincon')

% Solve the optimal control problem
Data = OCP_tempchange_fmincon(ip,'');

% Solve the ODEs with the optimal profiles 
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Tavg = cal_Tavg(T_fmincon(:,n1+1:end-1));
dTdt = cal_dTdt(t_fmincon,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_fmincon));
Data.e2 = cal_dTdt_2norm(t_fmincon,Tavg ,dTdt_obj);
disp(['e2_fmincon = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_fmincon,cal_Tavg(T_fmincon(:,n1+1:end-1)),ip)

end

%% Optimization using fmincon with a piecewise linear control via parallel computing
switch fmincon_par
case 'on'
disp('Solving the optimal control problem with fmincon via parallel computing')

% Solve the optimal control problem
Data = OCP_tempchange_fmincon_par(ip,'');

% Solve the ODEs with the optimal profiles 
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Tavg = cal_Tavg(T_fmincon(:,n1+1:end-1));
dTdt = cal_dTdt(t_fmincon,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_fmincon));
Data.e2 = cal_dTdt_2norm(t_fmincon,Tavg ,dTdt_obj);
disp(['e2_fmincon = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_fmincon,cal_Tavg(T_fmincon(:,n1+1:end-1)),ip)

end

%% Direct single shooting with CasADi
switch CasADi_ss
case 'on'
disp('Solving the optimal control problem with CasADi, single shooting')

% Solve the optimal control problem
Data = OCP_tempchange_CSD_ss(ip,'');

% Solve the ODEs with the optimal profiles 
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
Tavg = cal_Tavg(T_CSD(:,n1+1:end-1));
dTdt = cal_dTdt(t_CSD,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_CSD));
Data.e2 = cal_dTdt_2norm(t_CSD,Tavg ,dTdt_obj);
disp(['e2_CSD = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_CSD,S_CSD) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_CSD,cal_Tavg(T_CSD(:,n1+1:end-1)),ip)

end

%% Direct multiple shooting with CasADi
switch CasADi_ms
case 'on'
disp('Solving the optimal control problem with CasADi, multiple shooting')

% Solve the optimal control problem
Data = OCP_tempchange_CSD_ms(ip,'');

% Solve the ODEs with the optimal profiles 
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases(t,T,[Data.Tb_opt,Data.tb],ip), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
Tavg = cal_Tavg(T_CSD(:,n1+1:end-1));
dTdt = cal_dTdt(t_CSD,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_CSD));
Data.e2 = cal_dTdt_2norm(t_CSD,Tavg ,dTdt_obj);
disp(['e2_CSD = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_CSD,S_CSD) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_CSD,cal_Tavg(T_CSD(:,n1+1:end-1)),ip)

end

%% Simulation-based approach
switch DAE
case 'on'
disp('Solving the optimal control problem with the simulation-based technique')

% Load the Python file
py_file = 'DAE_TempChange.py';

% DAE Solver from Python
Tb0 = 0.1443;  % consistent initial condition
cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
output_py = pyrunfile(py_file,'output_MATLAB',x=dTdt_obj,y=Tb0);
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
Tavg = cal_Tavg(T_DAE(:,n1+1:end-1));
dTdt = cal_dTdt(t_DAE,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_DAE));
Data.e2 = cal_dTdt_2norm(t_DAE,Tavg ,dTdt_obj);
disp(['e2_DAE = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_DAE,S_DAE) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_DAE,cal_Tavg(T_DAE(:,n1+1:end-1)),ip)

end

%% Optimization using GEKKO with a piecewise linear control
switch GEKKO
case 'on'
disp('Solving the optimal control problem with GEKKO')

% Load the Python file
py_file = 'Gekko_TempChange.py';
Tb0 = [];
cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
output_py = pyrunfile(py_file,'output_MATLAB',x=dTdt_obj);
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
Tavg = cal_Tavg(T_Gekko(:,n1+1:end-1));
dTdt = cal_dTdt(t_Gekko,Tavg);
Data.RMSE = sqrt(sum((Tavg-T_target).^2)/length(t_Gekko));
Data.e2 = cal_dTdt_2norm(t_Gekko,Tavg ,dTdt_obj);
disp(['e2_GEKKO = ' , num2str(Data.e2)])

% Plot the result
figure; plot_interface(t_Gekko,S_Gekko) 
figure; plot_Tb(Data.tb,ip.temp_non(Data.Tb_opt))
figure; plot_dTdt(t_Gekko,cal_Tavg(T_Gekko(:,n1+1:end-1)),ip)

end