function xs = xs_generator(r,dim)
% This function defines the xs at a point in space defined by r

%This function is case dependent
%sigma = @(x,y) 5*x.^2 + 1;
sigma = @(x,y) 6 - 5*x.^2;
%sigma = @(x,y) 40*(x-0.5).^2 + 1;
%sigma = @(x,y) 10*sqrt(x) + 1;
%sigma = @(x,y) 10/(3-2*x);

xs = sigma(r(1),r(2));

