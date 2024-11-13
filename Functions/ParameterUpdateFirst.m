function [dn,dnu] = ParameterUpdateFirst(n0, nu0, w, Hxz0, Hy0)


    Hr=(nu0.')*(-w^(-2)*(Hxz0)+Hy0)*(nu0);
    gammar=( (1-eye(length(n0)) )./(n0-n0'+100*eye(length(n0)))).*Hr;
    
    dn=diag(Hr);
    dnu=-gammar*nu0;

end
