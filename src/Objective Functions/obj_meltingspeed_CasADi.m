function outputs = obj_meltingspeed_CasADi(t,S,dSdt_obj)

    import casadi.*
     
    dSdt = MX.zeros(length(t)-1,1);
    for i = 1:length(t)-1
        dt = t(i+1)-t(i);
        dSdt(i) = (S(i+1)-S(i))/dt;
    end
    
    outputs = sum((dSdt-dSdt_obj).^2);

end