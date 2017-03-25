function [ z ] = dvdlambda( r, NA, lambda, ni, z0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    k = 2*pi/lambda;
    v1 = @(rho) besselj(1, k * NA * r .* rho) .* sin(0.5 * k * rho.^2 * z0 * (NA^2 / ni)) .* (rho .^ 2);
    v2 = @(rho) besselj(0, k * NA * r .* rho) .* cos(0.5 * k * rho.^2 * z0 * (NA^2 / ni)) .* (rho .^ 3);
    temp1 = integral(v1,0,1,'ArrayValued',true,'AbsTol',1e-12,'RelTol',1e-9);
    temp2 = integral(v2,0,1,'ArrayValued',true,'AbsTol',1e-12,'RelTol',1e-9);
    z = k * (NA/lambda) * r .* temp1 - 0.5 * k * (z0/lambda) * (NA^2/ni) * temp2;

end

