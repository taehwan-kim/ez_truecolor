function p = make_disk(N,r)
%

x = zeros(N,N);
ncen = 1+N/2;

for i=1:N
    for j=1:N
        x(i,j) = sqrt((i-ncen)^2+(j-ncen)^2)/(N/2);
    end
end

p = x<=r;
