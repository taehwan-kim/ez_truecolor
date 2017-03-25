clc; clear; close all;

N = 1500;
NA = 1.4;
nm = 1.5;
z0 = 0e-9;
lambda = 600e-9;
pixelsize = 100e-9;
numofpixels = 10;
backgroundphotons = 10;
rep = 100;


pixelindex = linspace(0, (numofpixels/2)*pixelsize, numofpixels/2+1);

qraw = @(x, y) uz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2 + vz0(sqrt(x.^2+y.^2), NA, lambda, nm, z0).^2;

qarray = zeros(1, (numofpixels/2)^2);
psfmodel = zeros(numofpixels/2, numofpixels/2);
coordinatesx = zeros(numofpixels/2, numofpixels/2);
coordinatesy = zeros(numofpixels/2, numofpixels/2);
coordinates = zeros(numofpixels, numofpixels,2);

iter = ((numofpixels^2/4));

for i=1:iter

    k = mod(i, numofpixels/2);
    if (k==0) 
        k=numofpixels/2;
    end
    j = ((i-k)/(numofpixels/2))+1;

    qarray(i) = integral2(qraw,pixelindex(j),pixelindex(j+1),pixelindex(k),pixelindex(k+1));
    
end

Az0 = sum(qarray)*4;

muarray = qarray * N / Az0;

for i=1:numofpixels/2
    for j=1:numofpixels/2
        psfmodel(i,j) = muarray((i-1)*numofpixels/2 + j);
        coordinatesx(i,j) = pixelsize * j - pixelsize * 0.5;
        coordinatesy(i,j) = pixelsize * i - pixelsize * 0.5;
    end
end


psfmodel = [rot90(psfmodel,2) flipud(psfmodel) ; fliplr(psfmodel) psfmodel];
coordinatesx = [rot90(coordinatesx,2) flipud(coordinatesx) ; fliplr(coordinatesx) coordinatesx];
coordinatesy = [rot90(coordinatesy,2) flipud(coordinatesy) ; fliplr(coordinatesy) coordinatesy];
coordinates(:,:,1) = coordinatesx*1e9;
coordinates(:,:,2) = coordinatesy*1e9;

background = backgroundphotons * ones(size(psfmodel));

expected_size = 0.61e9*lambda/NA;

middlepatch = integral2(qraw,-1*pixelsize/2,1*pixelsize/2,-1*pixelsize/2,1*pixelsize/2)*N/Az0;

params = [middlepatch, 0, expected_size, 0, expected_size, backgroundphotons];


expectedsignal = zeros(numofpixels, numofpixels, rep);
estimated_size = zeros(1, rep);


for i=1:rep

    expectedsignal(:,:,i) = poissrnd(psfmodel) + poissrnd(background);
    solparz_pixel = lsqcurvefit(@d2gaussian,params,coordinates,expectedsignal(:,:,i));
    estimated_size(i) = sqrt(solparz_pixel(3)^2+solparz_pixel(5)^2);

end




