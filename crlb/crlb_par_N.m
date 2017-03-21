clc; close all;

NA = 1.25;
nm = 1.5;
lambda = 600e-9;
pixelsize = 100e-9;
numofpixels = 20;
background = 10;
z0 = 0;

N = 200*linspace(1,10,10);

pixelindex = linspace(0, (numofpixels/2)*pixelsize, numofpixels/2+1);

iter = ((numofpixels^2/4));

dmudxarray = zeros(length(N), iter);
dmudyarray = zeros(length(N), iter);
dmudlambdaarray = zeros(length(N), iter);
dmudNarray = zeros(length(N), iter);

I33 = zeros(3,3,length(N),iter);
sumI33 = zeros(3,3,length(N));
crlb = zeros(3,3,length(N));
sigmas = zeros(3,length(N));


for h=1:length(N)

    qraw = @(x, y) uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2 + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2;
    dAz0dlambdatemp = @(x, y) (2*uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).*dudlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0) + 2*vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).*dvdlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0));
    
    Az0 = integral2(qraw,-1*pixelindex(end),pixelindex(end),-1*pixelindex(end),pixelindex(end));
    dAz0dlambda = integral2(dAz0dlambdatemp,-1*pixelindex(end),pixelindex(end),-1*pixelindex(end),pixelindex(end));
    
    q = @(x, y) qraw(x, y) / Az0;
    
    dqdx = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dudx(x, y, NA, lambda, nm, z0) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dvdx(x, y, NA, lambda, nm, z0) ) ;
    dqdy = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dudx(y, x, NA, lambda, nm, z0) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dvdx(y, x, NA, lambda, nm, z0) ) ;
    dqdlambda = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dudlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0) .* dvdlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0) ) - (1 / Az0) * dAz0dlambda * q(x, y);

    parfor i=1:iter
    
        k = mod(i, numofpixels/2);
        if (k==0) 
            k=10;
        end
        j = ((i-k)/(numofpixels/2))+1;
    
        dmudxarray(h,i) = 2*N(h)*integral2(dqdx,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
        dmudyarray(h,i) = 2*N(h)*integral2(dqdy,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
        dmudlambdaarray(h,i) = 2*N(h)*integral2(dqdlambda,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
        dmudNarray(h,i) = integral2(q,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    
    end

    for i=1:iter
        
        I33(1,1,h,i) = (4/(background + N(h)*dmudNarray(h,i))) * dmudxarray(h,i) * dmudxarray(h,i);
        I33(2,2,h,i) = (4/(background + N(h)*dmudNarray(h,i))) * dmudyarray(h,i) * dmudyarray(h,i);
        I33(3,3,h,i) = (4/(background + N(h)*dmudNarray(h,i))) * dmudlambdaarray(h,i) * dmudlambdaarray(h,i);
        I33(1,2,h,i) = (4/(background + N(h)*dmudNarray(h,i))) * dmudxarray(h,i) * dmudyarray(h,i);
        I33(1,3,h,i) = (4/(background + N(h)*dmudNarray(h,i))) * dmudxarray(h,i) * dmudlambdaarray(h,i);
        I33(2,3,h,i) = (4/(background + N(h)*dmudNarray(h,i))) * dmudyarray(h,i) * dmudlambdaarray(h,i);
        I33(2,1,h,i) = I33(1,2,h,i);
        I33(3,1,h,i) = I33(1,3,h,i);
        I33(3,2,h,i) = I33(2,3,h,i);
        
    end
    for l=1:iter
        sumI33(:,:,h) = sumI33(:,:,h) + I33(:,:,h,l);
    end
    
    crlb(:,:,h) = sumI33(:,:,h)^-1;
    
    for l=1:3
        sigmas(l,h) = sqrt(crlb(l,l,h));
    end
end

dataforsave = zeros(3,length(N));
dataforsave(1,:) = N;
dataforsave(2,:) = sigmas(1,:);
dataforsave(3,:) = sigmas(3,:);

filename = ('par_N_data.mat');
save(filename, 'dataforsave');


%CRLB=I33^-1;

%sigmas = zeros(1,3);
%for i=1:3
    %sigmas(i) = sqrt(CRLB(i,i));
%end

