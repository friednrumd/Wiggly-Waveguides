function [dv, dn] = ParameterUpdateRK4(v0,n0,dH,w,dw,S0inv)

    [v1, n1]=ParameterUpdateFirst(v0              ,n0        ,dH(w)     );
    [v2, n2]=ParameterUpdateFirst(v0+S0inv*v1*dw/2,n0+n1*dw/2,dH(w+dw/2));
    [v3, n3]=ParameterUpdateFirst(v0+S0inv*v2*dw/2,n0+n2*dw/2,dH(w+dw/2));
    [v4, n4]=ParameterUpdateFirst(v0+S0inv*v3*dw  ,n0+n3*dw  ,dH(w+dw  ));


    dn=      1/6*(n4    +2*n3     +2*n2     +n1);
    dv=S0inv*1/6*(v4    +2*v3     +2*v2     +v1);

end

