function F = d2gaussian_zeromean(x,xdata)
 F = x(1)*exp(-((xdata(:,:,1)).^2/(2*x(2)^2) + (xdata(:,:,2)).^2/(2*x(3)^2) ));
 
 %% x: parameters
 % x1: size
 % x2: x axis stddev
 % x3: y axis stddev
 
 