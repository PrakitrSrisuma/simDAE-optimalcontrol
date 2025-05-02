function plot_dTdt(t,T,input)

    dTdt = zeros(length(t)-1,1);
    for i = 1:length(t)-1
        delt = t(i+1)-t(i);
        dTdt(i) = (T(i+1)-T(i))/delt;
    end
    plot(t(1:end-1),dTdt,'-b','linewidth',1.5) 
    ylabel('{\itd}\Theta_{avg}/{\itd}\tau')
    xlabel('\tau')
    ylim([input.dTdt_obj-0.5*input.dTdt_obj input.dTdt_obj+0.5*input.dTdt_obj])

return
