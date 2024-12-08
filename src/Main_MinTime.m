% ==============================================================================
% This is a top-level routine.
% Simulation-based approach for fast optimal control
% 
% This routine considers the optimal control problem that minimizes the thawing 
% time in cell thawing.
% ==============================================================================
close all; clear; clc;

%% Preparation
% Add paths
addpath('Saved Data','Python','Events','Input Data','Plotting', ...
    'Objective Functions','Optimal Control', 'PDEs','Calculations')
addpath('casadi-3.6.3-windows64-matlab2018b')  % add your path to CasADi here

% Typical input data
input_def = get_inputdata;  % unpack all default input data
input_def.nc = 1;  % control interval
input = input_processing(input_def);  % processing the input  
Tb_def = input.Tb;  % default heater temperature 
n1 = input.n1;
endtime = input.endtime_temp;  % final time
tspan = [(0:input.dt_plot:endtime)'; endtime];  % time span with interval = dt_plot
tspan = unique(tspan);  % remove duplicated value
T_ini = input.IC;  % initial condition

% ODE solvers
tol_ode = input.tol_ode;
option_ode1 = input.option_ode1;
option_ode2 = input.option_ode2;

% Optimal control
nrun = input.nrun;  % number of simulation runs
nc = input.nc;  % number of control intervals


%% Option selection ('on' or 'off')
orgsim = 'off';  % original simulation
Ipopt = 'off';  % Ipopt
fmincon = 'off';  % fmincon
fmincon_par = 'off';  % fmincon with parallel computing
CasADi_ss = 'off';  % CasADi single shooting
CasADi_ms = 'off';  % CasADi multiple shooting
simDAE = 'off';  % Simulation-based method using the GEKKO solver called from Python
GEKKO = 'off';  % optimization with GEKKO


%% ODE solver
switch orgsim
case 'on'
[t_sim,T_sim] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,0,input), tspan, T_ini, option_ode2);
S_sim = T_sim(:,end);  % interface position
dTdt = cal_dTdt(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)));

figure; plot_interface(t_sim,S_sim)
figure; plot_Tavg(t_sim,cal_Tavg(T_sim(:,n1+1:end-1)))

end


%% Opimization using Ipopt
switch Ipopt
case 'on'
disp('Solving the optimal control problem with Ipopt')

% Solve the optimal control problem
Data = OCP_mintime_Ipopt(input);

% Solve the ODEs with the optimal profiles
[t_Ipopt,T_Ipopt] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Data.Tb_opt(1),input), tspan, T_ini, option_ode2);
S_Ipopt = T_Ipopt(:,end);
Tavg = cal_Tavg(T_Ipopt(:,n1+1:end-1));
disp(['Tb_opt_Ipopt = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_Ipopt,S_Ipopt) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))

end


%% Opimization using fmincon with a piecewise linear control
switch fmincon
case 'on'
disp('Solving the optimal control problem with fmincon')

% Solve the optimal control problem
Data = OCP_mintime_fmincon(input);

% Solve the ODEs with the optimal profiles
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Data.Tb_opt(1),input), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Tavg = cal_Tavg(T_fmincon(:,n1+1:end-1));
disp(['Tb_opt_fmincon = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))


end


%% Opimization using fmincon with a piecewise linear control via parallel computing
switch fmincon_par
case 'on'
disp('Solving the optimal control problem with fmincon via parallel computing')

% Solve the optimal control problem
Data = OCP_mintime_fmincon_par(input);

% Solve the ODEs with the optimal profiles
[t_fmincon,T_fmincon] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Data.Tb_opt(1),input), tspan, T_ini, option_ode2);
S_fmincon = T_fmincon(:,end);
Tavg = cal_Tavg(T_fmincon(:,n1+1:end-1));
disp(['Tb_opt_fmincon = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_fmincon,S_fmincon) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))

end


%% Direct single shooting with CasADi
switch CasADi_ss
case 'on'
disp('Solving the optimal control problem with CasADi, single shooting')

% Solve the optimal control problem
Data = OCP_mintime_CSD_ss(input);

% Solve the ODEs with the optimal profiles
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Data.Tb_opt(1),input), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
Tavg = cal_Tavg(T_CSD(:,n1+1:end-1));
disp(['Tb_opt_CSDss = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_CSD,S_CSD) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))

end


%% Direct multiple shooting with CasADi
switch CasADi_ms
case 'on'
disp('Solving the optimal control problem with CasADi, multiple shooting')

% Solve the optimal control problem
Data = OCP_mintime_CSD_ms(input);

% Solve the ODEs with the optimal profiles
[t_CSD,T_CSD] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Data.Tb_opt(1),input), tspan, T_ini, option_ode2);
S_CSD = T_CSD(:,end);
Tavg = cal_Tavg(T_CSD(:,n1+1:end-1));
disp(['Tb_opt_CSDms = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_CSD,S_CSD) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))

end


%% Solving the optimal control problem with the simulation-based technique
switch simDAE
case 'on'
disp('Solving the optimal control problem with the simulation-based technique')

% Solving the DAE
tic; [t_DAE,T_DAE] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,0,input), tspan, T_ini, option_ode2); 
S_DAE = T_DAE(:,end);
Data.twall = toc;
Data.Tb_opt = [input.Tmax;input.Tmax];
Data.tb = [t_DAE(1); t_DAE(end)];
Data.ext = '';
Data.Tb_ini = '';
Data.RMSE = '';
Data.e2 = '';
disp(['Tb_opt_simDAE = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_DAE,S_DAE) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))

end

%% Solving the optimal control problem with  GEKKO optimization
switch GEKKO
case 'on'
disp('Solving the optimal control problem with GEKKO')

% Load the Python file
py_file = 'Gekko_MinTime.py';
Tb0 = [];
cd([fileparts(matlab.desktop.editor.getActiveFilename),'\Python'])
output_py = pyrunfile(py_file,'output_MATLAB');
Theta_opt_Gekko = double(output_py{2})';
Data.Tb_opt = input.temp_dim(Theta_opt_Gekko);
Data.tb = double(output_py{3})';
Data.twall = output_py{1};
Data.ext = '';
Data.Tb_ini = Tb0;
cd(fileparts(matlab.desktop.editor.getActiveFilename))
    
% Solve the ODEs with the optimal profiles from the simulation-based technique
[t_Gekko,T_Gekko] = ode15s(@(t,T)PDE_2Phases_Robin(t,T,Data.Tb_opt(1),input), tspan, T_ini, option_ode2);
S_Gekko = T_Gekko(:,end);
Data.RMSE = '';
Data.e2 = '';
disp(['Tb_opt_GEKKO = ' , num2str(Data.Tb_opt(1))])

% Plot the result
figure; plot_interface(t_Gekko,S_Gekko) 
figure; plot_Tb(Data.tb,input.temp_non(Data.Tb_opt))

end

