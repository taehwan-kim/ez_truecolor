clc; clear all; close all;

N = 100;
NA = 1.4;
nm = 1.5;
z0 = 0e-9;
lambda = 400e-9;
pixelsize = 100e-9;
numofpixels = 20;

pixelindex = linspace(0, (numofpixels/2)*pixelsize, numofpixels/2+1);

qraw = @(x, y) uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2 + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2;
dAz0dlambdatemp = @(x, y) (2*uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).*dudlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0) + 2*vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).*dvdlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0));

Az0 = integral2(qraw,-1*pixelindex(end),pixelindex(end),-1*pixelindex(end),pixelindex(end));
dAz0dlambda = integral2(dAz0dlambdatemp,-1*pixelindex(end),pixelindex(end),-1*pixelindex(end),pixelindex(end));

q = @(x, y) qraw(x, y) / Az0;

dqdx = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dudx(x, y, NA, lambda, nm, z0) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dvdx(x, y, NA, lambda, nm, z0) ) ;
dqdy = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dudx(y, x, NA, lambda, nm, z0) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dvdx(y, x, NA, lambda, nm, z0) ) ;
dqdlambda = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dudlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dvdlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0) ) - (1 / Az0) * dAz0dlambda * q(x, y);

dmudxarray = zeros(numofpixels/2, numofpixels/2);
dmudyarray = zeros(numofpixels/2, numofpixels/2);
dmudlambdaarray = zeros(numofpixels/2, numofpixels/2);
dmudNarray = zeros(numofpixels/2, numofpixels/2);

I33 = zeros(3,3);

iter = ((numofpixels^2/4));

for i=1:iter

    k = mod(i, numofpixels/2);
    if (k==0) 
        k=10;
    end
    j = ((i-k)/(numofpixels/2))+1;

    dmudxarray(j,k) = 2*N*integral2(dqdx,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    dmudyarray(j,k) = 2*N*integral2(dqdy,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    dmudlambdaarray(j,k) = 2*N*integral2(dqdlambda,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    dmudNarray(j,k) = integral2(q,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));

    I33(1,1) = I33(1,1) + (4/(N*dmudNarray(j,k))) * dmudxarray(j,k) * dmudxarray(j,k);
    I33(2,2) = I33(2,2) + (4/(N*dmudNarray(j,k))) * dmudyarray(j,k) * dmudyarray(j,k);
    I33(3,3) = I33(3,3) + (4/(N*dmudNarray(j,k))) * dmudlambdaarray(j,k) * dmudlambdaarray(j,k);
    I33(1,2) = I33(1,2) + (4/(N*dmudNarray(j,k))) * dmudxarray(j,k) * dmudyarray(j,k);
    I33(1,3) = I33(1,3) + (4/(N*dmudNarray(j,k))) * dmudxarray(j,k) * dmudlambdaarray(j,k);
    I33(2,3) = I33(2,3) + (4/(N*dmudNarray(j,k))) * dmudyarray(j,k) * dmudlambdaarray(j,k);

end

I33(2,1) = I33(1,2);
I33(3,1) = I33(1,3);
I33(3,2) = I33(2,3);

CRLB=I33^-1;

sigmas = zeros(1,3);
for i=1:3
    sigmas(i) = sqrt(CRLB(i,i));
end

