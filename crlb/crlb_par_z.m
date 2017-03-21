clc; close all;

N = 1500;
NA = 1.25;
nm = 1.5;
lambda = 600e-9;
pixelsize = 100e-9;
numofpixels = 20;
background = 10;

z0_range = 500e-9;
z0 = linspace(0, z0_range, 3);

pixelindex = linspace(0, (numofpixels/2)*pixelsize, numofpixels/2+1);

iter = ((numofpixels^2/4));

dmudxarray = zeros(length(z0), iter);
dmudyarray = zeros(length(z0), iter);
dmudlambdaarray = zeros(length(z0), iter);
dmudNarray = zeros(length(z0), iter);

I33 = zeros(3,3,length(z0),iter);
sumI33 = zeros(3,3,length(z0));
crlb = zeros(3,3,length(z0));
sigmas = zeros(3,length(z0));


for h=1:length(z0)

    qraw = @(x, y) uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)).^2 + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)).^2;
    dAz0dlambdatemp = @(x, y) (2*uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)).*dudlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) + 2*vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)).*dvdlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)));
    
    Az0 = integral2(qraw,-1*pixelindex(end),pixelindex(end),-1*pixelindex(end),pixelindex(end));
    dAz0dlambda = integral2(dAz0dlambdatemp,-1*pixelindex(end),pixelindex(end),-1*pixelindex(end),pixelindex(end));
    
    q = @(x, y) qraw(x, y) / Az0;
    
    dqdx = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) .* dudx(x, y, NA, lambda, nm, z0(h)) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) .* dvdx(x, y, NA, lambda, nm, z0(h)) ) ;
    dqdy = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) .* dudx(y, x, NA, lambda, nm, z0(h)) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) .* dvdx(y, x, NA, lambda, nm, z0(h)) ) ;
    dqdlambda = @(x, y) (2 / Az0) * ( uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) .* dudlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) .* dvdlambda(sqrt(x.^2+y.^2), NA, lambda, nm, z0(h)) ) - (1 / Az0) * dAz0dlambda * q(x, y);

    parfor i=1:iter
    
        k = mod(i, numofpixels/2);
        if (k==0) 
            k=10;
        end
        j = ((i-k)/(numofpixels/2))+1;
    
        dmudxarray(h,i) = 2*N*integral2(dqdx,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
        dmudyarray(h,i) = 2*N*integral2(dqdy,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
        dmudlambdaarray(h,i) = 2*N*integral2(dqdlambda,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
        dmudNarray(h,i) = integral2(q,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    
    end

    for i=1:iter
        
        I33(1,1,h,i) = (4/(background + N*dmudNarray(h,i))) * dmudxarray(h,i) * dmudxarray(h,i);
        I33(2,2,h,i) = (4/(background + N*dmudNarray(h,i))) * dmudyarray(h,i) * dmudyarray(h,i);
        I33(3,3,h,i) = (4/(background + N*dmudNarray(h,i))) * dmudlambdaarray(h,i) * dmudlambdaarray(h,i);
        I33(1,2,h,i) = (4/(background + N*dmudNarray(h,i))) * dmudxarray(h,i) * dmudyarray(h,i);
        I33(1,3,h,i) = (4/(background + N*dmudNarray(h,i))) * dmudxarray(h,i) * dmudlambdaarray(h,i);
        I33(2,3,h,i) = (4/(background + N*dmudNarray(h,i))) * dmudyarray(h,i) * dmudlambdaarray(h,i);
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




%CRLB=I33^-1;

%sigmas = zeros(1,3);
%for i=1:3
    %sigmas(i) = sqrt(CRLB(i,i));
%end

