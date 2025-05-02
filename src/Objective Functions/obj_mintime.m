function outputs = obj_mintime(Tb, tspan, T_ini, option, input)

    [t,T] = ode15s(@(t,T)PDE_2Phases(t,T,Tb,input), tspan, T_ini, option);
    outputs = t(end)^2;

end