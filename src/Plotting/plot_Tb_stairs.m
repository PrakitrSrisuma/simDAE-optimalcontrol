function plot_Tb_stairs(t,Tb)
    
    Tb = [Tb;Tb(end)];
    stairs(t,Tb,'-b','linewidth',2) 
    ylabel('Normalized heater temperature')
    xlabel('Dimensionless time')
    %title('Optimal heater temperature')

return
