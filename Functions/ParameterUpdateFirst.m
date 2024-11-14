function [dn,dnu] = ParameterUpdateFirst(n0, nu0, dH)

    Hr=(nu0.')*(dH)*(nu0);
    
    dn=diag(Hr);
    dnu=-(( (1-eye(length(n0)) )./(n0-n0.'+100*eye(length(n0)))).*Hr)*nu0;

end

