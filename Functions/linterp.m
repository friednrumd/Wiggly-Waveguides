function y = linterp(xarray,yarray,x)

xn=IndexFinder(xarray,x);

if xarray(xn)<x
    y=( yarray(xn)-yarray(xn-1) )/( xarray(xn)-xarray(xn-1) )*( x-xarray(xn-1) )+yarray(xn-1);
elseif xarray(xn)==x
    y=yarray(xn);
elseif xarray(xn)>x
    y=( yarray(xn+1)-yarray(xn) )/( xarray(xn+1)-xarray(xn) )*( x-xarray(xn) )+yarray(xn);

end

