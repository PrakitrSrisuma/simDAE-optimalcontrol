function plot_abserror(t,S,S_target)

    error = abs(S-S_target);
    plot(t,error,'-b','linewidth',2) 
    ylabel('Normalized interface position')
    xlabel('Dimensionless time')
    %title('Optimal interface position')

return
