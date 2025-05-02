# Operating Conditions and General Parameters
b = 5e-3;  # radius (m)
q = 5800;  # heat flux (W)
Tini = -80;  # initial temperature (degC)
T0 = 37;  # boundary temperature (degC)
Tm = -2;  # melting point
TAH = 90;  # temperature alarm high
U = 850;  # heat transfer coefficient

# Fluid Properties
k1 = 2.22;  # thermal conductivity of ice (W/mK)
k2 = 0.556;  # thermal conductivity of water (W/mK)
Cp1 = 2.050;  # heat capacity of ice (kJ/kgK)
Cp2 = 4.2;  # heat capacity of water (kJ/kgK)
rho1 = 916;  # density of ice (kg/m3)
rho2 = 1000;  # density of water (kg/m3)
L = 334;  # Latent heat of melting (kJ/kg)
alp1 = k1/(Cp1*10e3*rho1);  # thermal diffusivity of ice (m2/s)
alp2 = k2/(Cp2*10e3*rho2);   # thermal diffusivity of water (m2/s)
rho3 = 0.5*(rho1+rho2);
alp = alp2/alp1;  # Dimensionless thermal diffusivity

## Dimensionless Parameters
tempconvert = T0-Tm;  # temperature conversion
timeconvert = alp1*60/(b)**2;  # time conversion (from actual to dimensionless)
T0_d = 1;  # dimensionless boundary temperature
Tm_d = 0;  # dimensionless melting point
Ste = Cp1*tempconvert/L;  # Stefan number (dimensionless)
Tini_d = (Tini-Tm)/tempconvert;  # dimensionless initial temperature
TAH_d = (TAH-Tm)/tempconvert;  # dimensionless temperature alarm high
U_d = U*b/k2;

## Discretization and Geometry
n1 = 20;  # r-direction
n2 = 30;  # theta-direction
endtime = 7;
dR = 1/n1;

## Control
dSdt_obj = -0.1;
dTdt_obj = 0.04;

#Others
tol = 1e-4;  # Stefan boundary tolerance (to avoid singularity)
