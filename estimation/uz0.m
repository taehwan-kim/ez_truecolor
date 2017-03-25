function [ y ] = uz0( r, NA, lambda, ni, z0 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    k = 2*pi/lambda;
    u = @(rho) besselj(0, k * (NA) * r .* rho) .* cos(0.5 * k * rho.^2 * z0 * (NA^2/ni) ) .* rho;
    y = integral(u,0,1,'ArrayValued',true,'AbsTol',1e-12,'RelTol',1e-9);

end

