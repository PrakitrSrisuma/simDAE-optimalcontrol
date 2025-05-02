function plot_dSdt(t,S,input)

    dSdt = zeros(length(t)-1,1);
    for i = 1:length(t)-1
        delt = t(i+1)-t(i);
        dSdt(i) = (S(i+1)-S(i))/delt;
    end
    plot(t(1:end-1),dSdt,'-b','linewidth',1.5) 
    ylabel('{\itdS}/{\itd}\tau')
    xlabel('\tau')
    ylim([input.dSdt_obj+0.5*input.dSdt_obj input.dSdt_obj-0.5*input.dSdt_obj])
    %title('Optimal interface position')

return
