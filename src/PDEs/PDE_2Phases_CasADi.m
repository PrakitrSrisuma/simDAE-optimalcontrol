function outputs = PDE_2Phases_CasADi(T,t,t_end,Tb,input)
import casadi.*

n1 = input.n1;
dR = input.dR;
Ste = input.Ste;
Tm = input.Tm_d;
alp = input.alp;
k1 = input.k1;
k2 = input.k2;
U = input.U_d;
Tb = input.temp_non(t*Tb(2)/t_end + (t_end-t)*Tb(1)/t_end);

%% Model
%\ Prepare vectors
diagonal1 = MX.zeros(n1,1);
subdiag1 = MX.zeros(n1,1);
supdiag1 = MX.zeros(n1,1);
diagonal2 = MX.zeros(n1,1);
subdiag2 = MX.zeros(n1,1);
supdiag2 = MX.zeros(n1,1);
RHS1 = MX.zeros(n1,1);
RHS2 = MX.zeros(n1,1);

%\ Stefan Boundary
S = T(end);  % Stefan boundary
dSdt = (Ste/(S))*(1.5*0-2*T(n1)/dR+T(n1-1)/(2*dR))-((k2*Ste)/(k1*(1-S)))*(-1.5*0+2*T(n1+1)/dR-T(n1+2)/(2*dR));
%dSdt = min(dSdt,0);

% S = if_else(S <= input.tol,input.tol,S);

%\ Solid Phase
% Prepare variables
for i = 1:n1
    a1 = 1/((S*dR)^2);  % coefficient a1
    b1 = 1/(2*i*(S*dR)^2);  % coefficient b1  
    c1 = 1/(2*S);  % coefficient c1
    
	diagonal1(i,1) = -2*a1;  % diagonal elements
	subdiag1(i,1) = a1-b1-c1*dSdt;  % subdiagonal elements
	supdiag1(i,1) = a1+b1+c1*dSdt;  % superdiagonal elements
    
    % Boundary condition
    if i == n1-1
        RHS1(n1,1) = (a1+b1+c1*dSdt)*(Tm);  % at the Stefan boundary
    end
end

diagonal1(1) = -4*a1;  % diagonal element handling the symmetry condition
supdiag1 = [4*a1; supdiag1(1:end-2)];  % superdiagonal element vector
subdiag1 = subdiag1(1:end-1);  % subdiagonal element vectpr

% Construct a system of equations
A1 = diag(diagonal1);

for i = 1:n1-1
    A1(i,i+1) = supdiag1(i);
end

for i = 2:n1
    A1(i,i-1) = subdiag1(i-1);
end

dTdA1 = A1*(T(1:n1))+RHS1;


%\ Liquid Phase
% Prepare variables
for j = 1:n1
    a2 = alp/(((1-S)^2)*(dR^2));  % coefficient a2
    b2 = alp/(2*(1-S)*(S+j*dR*(1-S))*dR);  % coefficient b2 
    c2 = (1-j*dR)/((1-S)*2*dR);  % coefficient c2
    
	diagonal2(j,1) = -2*a2;  % diagonal elements
	subdiag2(j,1) = a2-b2-c2*dSdt;  % subdiagonal elements
	supdiag2(j,1) = a2+b2+c2*dSdt;  % superdiagonal elements
    
    % Boundary conditions
    if j == 1
        RHS2(1,1) = (a2-b2-c2*dSdt)*(Tm);  % at the Stefan boundary
    end
    
    if j == n1
        d2 = 2*dR*U*(1-S);
        diagonal2(j) = diagonal2(j) - (a2+b2+c2*dSdt)*d2;
        subdiag2(j) = subdiag2(j) + (a2+b2+c2*dSdt);
        RHS2(n1,1) = (a2+b2+c2*dSdt)*d2*Tb;  % at the outer boundary
    end
    
end

subdiag2 = subdiag2(2:end);  % subdiagonal element vector
supdiag2 = supdiag2(1:end-1);  % superdiagonal element vector

% Construct a system of equations
A2 = diag(diagonal2);

for i = 1:n1-1
    A2(i,i+1) = supdiag2(i);
end

for i = 2:n1
    A2(i,i-1) = subdiag2(i-1);
end

dTdA2 = A2*(T(n1+1:2*n1)) + RHS2;

% dTdA1 = if_else(S <= input.tol,0,dTdA1);
% dTdA2 = if_else(S <= input.tol,0,dTdA2);
% dSdt = if_else(S <= input.tol,0,dSdt);

% Model dynamic
dTdt_model = [dTdA1; dTdA2; dSdt];

%% Export
outputs = dTdt_model;

return