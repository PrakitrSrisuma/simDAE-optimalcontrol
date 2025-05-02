function outputs = cal_dTdt_2norm(t,T,dTdt_obj)

nt = length(t)-1;
dTdt = zeros(nt,1);
for i = 1:length(t)-1
    delt = t(i+1)-t(i);
    dTdt(i) = (T(i+1)-T(i))/delt;
end

e2 = dTdt-dTdt_obj*ones(nt,1);
outputs = (1/sqrt(nt))*norm(e2,2);

return