function xs = xs_generator(r,dim)
% This function defines the xs at a point in space defined by r

%This function is case dependent
xs1 = 1.0;
xs2 = 10.0;

%Bottom and top middle
%if(r(1)/dim(1) > 1/3 && r(1)/dim(1) <= 2/3) && (r(2)/dim(2) > 1/3 && r(2)/dim(2) <= 2/3)
if(r(1)/dim(1) > 1/4 && r(1)/dim(1) <= 3/4) && (r(2)/dim(2) > 1/4 && r(2)/dim(2) <= 3/4)
%if(r(1)/dim(1) + (2/20)/dim(1) > 1/4 && r(1)/dim(1) + (2/20)/dim(1) <= 3/4) && (r(2)/dim(2) > 1/4 && r(2)/dim(2) <= 3/4)    %trying to offset the pin 
    xs = xs2;
else
    xs = xs1;
end


