function outputs = get_inputdata
 
% Optimal Control-related Parameters
input.nc = 1;  % number of control intervals
input.nrun = 1;  % number of simulation runs
input.dSdt_obj = -0.1;  % target melting speed
input.dTdt_obj = 0.04;  % target temperature change
input.dT_obj = 0.1;  % target driving force

% Operating Conditions and General Parameters
input.b = 5e-3;  % radius (m)
input.q = 5800;  % heat flux (W)
input.Tini = -80;  % initial temperature (degC)
input.T0 = 37;  % reference temperature (degC)
input.Tm = -2;  % melting point
input.TAH = 90;  % temperature alarm high
input.U = 850;  % heat transfer coefficient
input.Tb = 37;  % heater temperature (degC)
input.alp = [];

% Fluid Properties
input.k1 = 2.22;  % thermal conductivity of ice (W/mK)
input.k2 = 0.556;  % thermal conductivity of water (W/mK)  
input.Cp1 = 2.050;  % heat capacity of ice (kJ/kgK)
input.Cp2 = 4.2;  % heat capacity of water (kJ/kgK)
input.rho1 = 916;  % density of ice (kg/m3) 
input.rho2 = 1000;  % density of water (kg/m3) 
input.L = 334;  % Latent heat of melting (kJ/kg)

% Numerics
input.n1 = 20;  % r-direction discretization
input.endtime = 1e5; % time step (min)
input.endtime_Itf = 7;
input.endtime_temp = 7;
input.tol = 1e-4;  % Stefan boundary tolerance (to avoid singularity)
input.tol_ode = 1e-10;  % tolerance for ODE solvers
input.tol_opt = 1e-7;  % tolerance for optimizer
input.n_it = 2000;  % number of iterations for optimizer
input.dt_plot = .05;  % for plotting and RMSE check
 
outputs = input;

return
