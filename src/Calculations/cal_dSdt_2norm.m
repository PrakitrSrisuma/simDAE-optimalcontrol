function outputs = cal_dSdt_2norm(t,S,dSdt_obj)

nt = length(t)-1;
dSdt = zeros(nt,1);
for i = 1:length(t)-1
    delt = t(i+1)-t(i);
    dSdt(i) = (S(i+1)-S(i))/delt;
end

e2 = dSdt-dSdt_obj*ones(nt,1);
outputs = (1/sqrt(nt))*norm(e2,2);

return