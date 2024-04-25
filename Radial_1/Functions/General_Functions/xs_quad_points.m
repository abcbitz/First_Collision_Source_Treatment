function xs = xs_quad_points(zn,r,dim)
% This function defines the xs at the quadrature points inside the zone. If
% actually doing the ray trace, these xs are just the zones xs. If
% integrating the actual function, the xs will differ across the zone, and
% this is needed.

global use_func;
if use_func == 1
    xs = xs_generator(r,dim);
else
    xs = zn.sigma_a;
end


