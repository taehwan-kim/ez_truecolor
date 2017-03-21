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

for i=1:numofpixels/2
    for j=1:numofpixels/2
        dmudxarray(i,j) = 2*N*integral2(dqdx,pixelindex(i),pixelindex(i+1),pixelindex(j),pixelindex(j+1));
        dmudyarray(i,j) = 2*N*integral2(dqdy,pixelindex(i),pixelindex(i+1),pixelindex(j),pixelindex(j+1));
        dmudlambdaarray(i,j) = 2*N*integral2(dqdlambda,pixelindex(i),pixelindex(i+1),pixelindex(j),pixelindex(j+1));
        dmudNarray(i,j) = integral2(q,pixelindex(i),pixelindex(i+1),pixelindex(j),pixelindex(j+1));
        
        I33(1,1) = I33(1,1) + (4/(N*dmudNarray(i,j))) * dmudxarray(i,j) * dmudxarray(i,j);
        I33(2,2) = I33(2,2) + (4/(N*dmudNarray(i,j))) * dmudyarray(i,j) * dmudyarray(i,j);
        I33(3,3) = I33(3,3) + (4/(N*dmudNarray(i,j))) * dmudlambdaarray(i,j) * dmudlambdaarray(i,j);
        I33(1,2) = I33(1,2) + (4/(N*dmudNarray(i,j))) * dmudxarray(i,j) * dmudyarray(i,j);
        I33(1,3) = I33(1,3) + (4/(N*dmudNarray(i,j))) * dmudxarray(i,j) * dmudlambdaarray(i,j);
        I33(2,3) = I33(2,3) + (4/(N*dmudNarray(i,j))) * dmudyarray(i,j) * dmudlambdaarray(i,j);
%         
    end
end

%  imagesc(dmudNarray);
A = dmudNarray/dmudNarray(1,1);
B = fliplr(A);
C = flipud(B);
D = flipud(A);
aa = [C D;B A];
imagesc(aa);
% 
% psf = PSFGenerator.get;
% temptemp=psf(:,:,14);

I33(2,1) = I33(1,2);
I33(3,1) = I33(1,3);
I33(3,2) = I33(2,3);
% 
CRLB=I33^-1;
% xacc = CRLB(1,1);
sigmas = zeros(1,3);
for i=1:3
    sigmas(i) = sqrt(CRLB(i,i));
end


% 
% muk = dmudNarray * N;



% temp = integral2(dqdx,0,pixelsize,0,pixelsize);