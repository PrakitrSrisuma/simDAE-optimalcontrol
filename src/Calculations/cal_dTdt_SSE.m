function outputs = cal_dTdt_SSE(t,T,dTdt_obj)

nt = length(t)-1;
dTdt = zeros(nt,1);
for i = 1:length(t)-1
    delt = t(i+1)-t(i);
    dTdt(i) = (T(i+1)-T(i))/delt;
end

outputs = sum((dTdt-dTdt_obj).^2);

return