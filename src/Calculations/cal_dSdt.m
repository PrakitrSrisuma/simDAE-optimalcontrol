function outputs = cal_dSdt(t,S)

nt = length(t)-1;
dSdt = zeros(nt,1);
for i = 1:length(t)-1
    delt = t(i+1)-t(i);
    dSdt(i) = (S(i+1)-S(i))/delt;
end

outputs = dSdt;

return