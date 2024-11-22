function [dv,dn] = ParameterUpdateFirst(v0, n0, dHr)

    Hr=(v0.')*(dHr)*(v0);
    
    dn=diag(Hr);
    dv=-(( (1-eye( length(n0) ) )./(n0-n0.'+10*eye(length(n0)))).*Hr)*v0;

end

