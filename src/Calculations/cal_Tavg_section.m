function outputs = cal_Tavg_section(T,input)

Tavg = zeros(height(T),1);
R = (0+input.dR:input.dR:1);
dA = sum(R);
for i = 1:length(Tavg)
    Tavg(i) = sum(R.*T(i,:))/dA;
end

outputs = Tavg;

return