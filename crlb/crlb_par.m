clc; clear all; close all;

N = 1500;
NA = 1.25;
nm = 1.5;
z0 = 0e-9;
lambda = 600e-9;
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

dmudxarray = zeros(1, (numofpixels/2)^2);
dmudyarray = zeros(1, (numofpixels/2)^2);
dmudlambdaarray = zeros(1, (numofpixels/2)^2);
dmudNarray = zeros(1, (numofpixels/2)^2);



iter = ((numofpixels^2/4));

I33 = zeros(3,3,iter);

parfor i=1:iter

    k = mod(i, numofpixels/2);
    if (k==0) 
        k=10;
    end
    j = ((i-k)/(numofpixels/2))+1;

    dmudxarray(i) = 2*N*integral2(dqdx,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    dmudyarray(i) = 2*N*integral2(dqdy,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    dmudlambdaarray(i) = 2*N*integral2(dqdlambda,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    dmudNarray(i) = integral2(q,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));

end

for i=1:iter
    
    I33(1,1,i) = (4/(N*dmudNarray(i))) * dmudxarray(i) * dmudxarray(i);
    I33(2,2,i) = (4/(N*dmudNarray(i))) * dmudyarray(i) * dmudyarray(i);
    I33(3,3,i) = (4/(N*dmudNarray(i))) * dmudlambdaarray(i) * dmudlambdaarray(i);
    I33(1,2,i) = (4/(N*dmudNarray(i))) * dmudxarray(i) * dmudyarray(i);
    I33(1,3,i) = (4/(N*dmudNarray(i))) * dmudxarray(i) * dmudlambdaarray(i);
    I33(2,3,i) = (4/(N*dmudNarray(i))) * dmudyarray(i) * dmudlambdaarray(i);
    I33(2,1,i) = I33(1,2,i);
    I33(3,1,i) = I33(1,3,i);
    I33(3,2,i) = I33(2,3,i);
    
end


%CRLB=I33^-1;

%sigmas = zeros(1,3);
%for i=1:3
    %sigmas(i) = sqrt(CRLB(i,i));
%end

