close all;
clear all;

% assume that NA = 1 and f-number = 0.5;
% 0.61 lambda = first dark ring = 2.5*FWHM
% wavelength: orange = 573, red = 702
% disk size: orange = 48, red = 28

%% input
lambda = 0.450;
disksize = 0.048;
pixelsize = 0.16;
N = 1024;                      % resolution: 2/N

%% computation of the 1D airy pattern
firstdark = lambda * 0.61;
FWHM = firstdark / 2.5;
k = N * FWHM;

pupil = make_pupil(N,2*k);    % construct the pupil function
rmax = N/(4*k);
w = zeros(N,N);               % perfect wavefront
z = psf(pupil,w);             % generate the spread function
ncen = 1+N/2;

% disk = make_disk(N,disksize);
% disk = disk .* exp(2*pi*j*w);

% zp = z(ncen:ncen+N/2-1,ncen);    % take the central portion
% xp = (0:N/2-1)/(N/2);            % generate x-coordinates

% calculate exact airy function
xe = linspace(0,2*N/2/(2*k),N);
ze = somb(2*xe);
ze = ze.*ze;

% calculate bead pattern modeled as a simple disk
disksize = disksize/2;
dsk = zeros(size(xe));
xe2 = (0:2/(length(xe)-1):2);
for i=1:N
    if (xe2(i)<=disksize)
        dsk(i) = 1;
    end
end

% calculate the pixel size
blocksize = 2 * floor(pixelsize*256) + 1;

%% rotate to create 2D pattern
zes = [flip(ze) ze(2:end)];
dsks = [flip(dsk) dsk(2:end)];
% V = [ zeros(1,20) ones(1,15) zeros(1,20)]; % Any symmetrical vector

V = zes;
Vdsk = dsks;
n = floor(numel(V)/2); 

r = [n:-1:0, 1:n]; % A vector of distance (measured in pixels) from the center of vector V to each element of V

% Now find the distance of each element of a square 2D matrix from it's centre. @(x,y)(sqrt(x.^2+y.^2)) is just the Euclidean distance function. 
ri = bsxfun(@(x,y)(sqrt(x.^2+y.^2)),r,r');

% Now use those distance matrices to interpole V
img = interp1(r(1:n+1),V(1:n+1),ri);
imgdsk = interp1(r(1:n+1),Vdsk(1:n+1),ri);

% The corners will contain NaN because they are further than any point we had data for so we get rid of the NaNs
img(isnan(img)) = 0; % or instead of zero, whatever you want your background colour to be
imgdsk(isnan(imgdsk)) = 0;

temp = size(img);
imgcen = (temp(1)+1)/2;
img = img(imgcen-imgcen/2:imgcen+imgcen/2, imgcen-imgcen/2:imgcen+imgcen/2);
imgdsk = imgdsk(imgcen-imgcen/2:imgcen+imgcen/2, imgcen-imgcen/2:imgcen+imgcen/2);


%% convolve psf with the disk image
output = conv2(imgdsk, img, 'same');
output = 20*output/max(max(output)); % normalize the convolution

%% pixelation
fun = @(block_struct) ...
   mean2(block_struct.data) * ones(size(block_struct.data));
imgpixel = blockproc(img,[blocksize, blocksize],fun);
outputpixel = blockproc(output,[blocksize, blocksize],fun);


%% plot the results
% figure(1);
% plot(xe2,ze,xe2,dsk);
% xlabel('radius');
% ylabel('irradiance');
% 
% figure(2);
% x = (-1:2/(N-1):1);
% y = (-1:2/(N-1):1);
% imagesc(x,y,log_image(img,5));

% figure(3);
% x = (-1:2/(N-1):1);
% y = (-1:2/(N-1):1);
% imagesc(x,y,output);

% generate and save the images.
% Use logrithmic transformation for output image
% figure(3);
% imshow(log_image(imgdsk,5));
% figure(4);
% imshow(log_image(img,5));
% figure(5);
% imshow(log_image(output,5));
% imwrite(255*log_image(img,5),gray(255),'airy.bmp');
% imwrite(255*pupil,gray(255),'pupil.bmp');

%% gaussian fit of the airy disk

params = [1, 0, 0.5, 0, 0.5];
% params_zeromean = [1, 0.5, 0.5];

 % parameters
 % x1: size
 % x2: x axis mean
 % x3: x axis stddev
 % x4: y axis mean
 % x5: y axis stddev

[X, Y] = meshgrid(-1:2/(N):1);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;

solparz = lsqcurvefit(@d2gaussian,params,xdata,img);
solparz_pixel = lsqcurvefit(@d2gaussian,params,xdata,imgpixel);

% solparout = lsqcurvefit(@d2gaussian,params,xdata,output);
% solparout_pixel = lsqcurvefit(@d2gaussian,params,xdata,outputpixel);

% solparout_zeromean = lsqcurvefit(@d2gaussian_zeromean,params_zeromean,xdata,output);

solz = solparz(3);
solz_pixel = solparz_pixel(3);
% sol = solparout(3);
% sol_pixel = solparout_pixel(3);
% sol_zeromean = solparout_zeromean(3);
