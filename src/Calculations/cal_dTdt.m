function outputs = cal_dTdt(t,T)

nt = length(t)-1;
dTdt = zeros(nt,1);

for i = 1:nt
    delt = t(i+1)-t(i);
    dTdt(i) = (T(i+1)-T(i))/delt;
end

outputs = dTdt;

return