function outputs = obj_meltingspeed_pc(Tb, tb, dSdt_obj, T_ini, option, input)

    tf = [];
    S = [];

    for i = 1:input.nc
        tspan = (tb(i):input.dt_plot:tb(i+1));
        [t,T] = ode15s(@(t,T)PDE_2Phases(t,T,Tb(i),input), tspan, T_ini, option);
        if i~=input.nc
            tf = [tf;t(1:end-1)];
            S = [S;T(1:end-1,end)];
        else
            tf = [tf;t];
            S = [S;T(:,end)];
        end
        T_ini = T(end,:);
    end

    outputs = cal_dSdt_SSE(tf,S,dSdt_obj);

end