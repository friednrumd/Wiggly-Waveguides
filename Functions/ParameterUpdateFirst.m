function [dn,dnu] = ParameterUpdateFirst(n0, nu0, w, Hxz0, Hy0)


    Hr=(nu0.')*(-(Hxz0)*w^(-2)+Hy0)*(nu0);
    gammar=((1-eye(length(n0)))./(n0-n0'+10^-20)).*Hr;
    
    dn=diag(Hr);
    dnu=-gammar*nu0;

end

