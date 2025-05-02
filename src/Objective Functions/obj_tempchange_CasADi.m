function outputs = obj_tempchange_CasADi(t,T,dTdt_obj)

import casadi.*
 
nt = length(t);
dTdt = MX.zeros(nt-1,1);
for i = 1:nt-1
    delt = t(i+1)-t(i);
    dTdt(i) = (T(i+1)-T(i))/delt;
end

outputs = sum((dTdt-dTdt_obj).^2);


end