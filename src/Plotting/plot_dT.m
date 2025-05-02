function plot_dT(t,T,Tb,input)

    plot(t,T,'-b','linewidth',2) 
    hold on 
    plot(t,Tb,'--r','linewidth',2)
    ylabel('Normalized temperature')
    xlabel('Dimensionless time')
    legend('Average temperature','Heater temperature')
    % ylim([input.dTdt_obj-0.5*input.dTdt_obj input.dTdt_obj+0.5*input.dTdt_obj])

return
