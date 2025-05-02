function outputs = cal_dSdt_SSE(t,S,dSdt_obj)

nt = length(t)-1;
dSdt = zeros(nt,1);
for i = 1:length(t)-1
    delt = t(i+1)-t(i);
    dSdt(i) = (S(i+1)-S(i))/delt;
end

outputs = sum((dSdt-dSdt_obj).^2);

return