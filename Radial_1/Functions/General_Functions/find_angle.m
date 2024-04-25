function angle = find_angle(r,r_s)

x = r(1) - r_s(1);
y = r(2) - r_s(2);

if x >= 0.0
    angle = atan(y/x);  %quads 1 and 4
else
    if y >= 0.0
        angle = acos(x/sqrt(x^2+y^2));  %quad 2
    else
        angle = atan(y/x) + pi; %quad 3
    end
end

if angle < 0.0
    angle = 2*pi + angle;
end