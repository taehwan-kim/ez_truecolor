function [ z ] = dvdx( x, y, NA, lambda, ni, z0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    k = 2*pi/lambda;
    r = sqrt(x.^2 + y.^2);
    v = @(rho) besselj(1, k * NA * r .* rho) .* sin(0.5 * k * rho.^2 * z0 * (NA^2/ni)) .* (rho .^ 2);
    z = k * NA * x .* integral(v,0,1,'ArrayValued',true,'AbsTol',1e-12,'RelTol',1e-9) ./ r;

end