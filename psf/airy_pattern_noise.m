close all;
clear all;

% assume that NA = 1 and f-number = 0.5
% 0.61 lambda = first dark ring = 2.5*FWHM
% wavelength: orange = 573, red = 702
% disk size: orange = 48, red = 28

%% input
lambda = 0.573;
pixelsize = 0.16;
N = 1024;                      % resolution: 2/N
noise = 0.0;
gain = 2;

%% computation of the 1D airy pattern
firstdark = lambda * 0.61;
FWHM = firstdark / 2.5;
k = N * FWHM;

% disk = make_disk(N,disksize);
% disk = disk .* exp(2*pi*j*w);

% zp = z(ncen:ncen+N/2-1,ncen);    % take the central portion
% xp = (0:N/2-1)/(N/2);            % generate x-coordinates

% calculate exact airy function
xe = linspace(0,2*N/2/(2*k),N);
ze = somb(2*xe);
ze = ze.*ze;

% calculate bead pattern modeled as a simple disk

% calculate the pixel size
blocksize = 2 * floor(pixelsize*256) + 1;

%% rotate to create 2D pattern
zes = [flip(ze) ze(2:end)];
% V = [ zeros(1,20) ones(1,15) zeros(1,20)]; % Any symmetrical vector

V = zes;
n = floor(numel(V)/2); 

r = [n:-1:0, 1:n]; % A vector of distance (measured in pixels) from the center of vector V to each element of V

% Now find the distance of each element of a square 2D matrix from it's centre. @(x,y)(sqrt(x.^2+y.^2)) is just the Euclidean distance function. 
ri = bsxfun(@(x,y)(sqrt(x.^2+y.^2)),r,r');

% Now use those distance matrices to interpole V
img = interp1(r(1:n+1),V(1:n+1),ri);

% The corners will contain NaN because they are further than any point we had data for so we get rid of the NaNs
img(isnan(img)) = 0; % or instead of zero, whatever you want your background colour to be

temp = size(img);
imgcen = (temp(1)+1)/2;
img = img(imgcen-imgcen/2:imgcen+imgcen/2, imgcen-imgcen/2:imgcen+imgcen/2);

%% monte carlo start

rep = 1;
width = zeros(1,rep);
inten = zeros(1,rep);

for i=1:rep

    %% add noise
    [X, Y] = meshgrid(-1:2/(N):1);
    imgtemp=gain*img+noise*(rand(size(X,1),size(Y,2)));

    %% pixelation
    fun = @(block_struct) ...
       mean2(block_struct.data) * ones(size(block_struct.data));
    imgpixel = blockproc(imgtemp,[blocksize, blocksize],fun);

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

    % solparz = lsqcurvefit(@d2gaussian,params,xdata,img);
    solparz_pixel = lsqcurvefit(@d2gaussian,params,xdata,imgpixel);


    % solparout_zeromean = lsqcurvefit(@d2gaussian_zeromean,params_zeromean,xdata,output);

    % solz = solparz(3);
    width(i) = sqrt(solparz_pixel(3)^2 + solparz_pixel(5)^2);
    inten(i) = solparz_pixel(1);
    % sol = solparout(3);
    % sol_pixel = solparout_pixel(3);
    % sol_zeromean = solparout_zeromean(3);
    
end

% hist(width);
% 
% savefile_width = './width_03_10.mat';
% savefile_inten = './inten_03_10.mat';
% save(savefile_width, 'width');
% save(savefile_inten, 'inten');
