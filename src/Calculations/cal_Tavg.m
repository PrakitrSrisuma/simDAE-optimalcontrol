function outputs = cal_Tavg(T)

Tavg = zeros(height(T),1);
for i = 1:length(Tavg)
    Tavg(i) = mean(T(i,:));
end

outputs = Tavg;

return