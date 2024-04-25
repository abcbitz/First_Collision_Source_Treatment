function [r,w] = quad_points_surf(grid,i,j,face,N,quad_method)

dx = grid(i,j).edge_lengths(1);
dy = grid(i,j).edge_lengths(2);

% digitsOld = digits(50);
% left = double(grid(i,j).location(1) - vpa(dx/2));
% right = double(grid(i,j).location(1) + vpa(dx/2));
% bottom = double(grid(i,j).location(2) - vpa(dy/2));
% top = double(grid(i,j).location(2) + vpa(dy/2));
% digits(digitsOld)

left = grid(i,j).location(1) - dx/2;
right = grid(i,j).location(1) + dx/2;
bottom = grid(i,j).location(2) - dy/2;
top = grid(i,j).location(2) + dy/2;

w = 0.0;

if strcmp(quad_method,'Gauss') == 1
    [y,w_v] = lgwt(N,bottom,top);
    [x,w_h] = lgwt(N,left,right);
    
    if face == 1
        r = [left*ones(N,1), flip(y)];
        w = w_v;
    elseif face == 2
        r = [flip(x), bottom*ones(N,1)];
        w = w_h;
    elseif face == 3
        r = [right*ones(N,1), flip(y)];
        w = w_v;
    elseif face == 4
        r = [flip(x), top*ones(N,1)];
        w = w_h;
    end
    
elseif strcmp(quad_method,'Gauss-Lobatto') == 1
    [y,w_v,~] = lglnodes(N-1); 
    w_v = w_v*dy/2; %dy/2 is the rescale factor from 2 (-1,1) to dy
    [x,w_h] = lglnodes(N-1);
    w_h = w_h*dx/2;
    
    if face == 1
        r = [left*ones(N,1), flip(y)*dy/2 + grid(i,j).location(2)];
        w = w_v;
    elseif face == 2
        r = [flip(x)*dx/2 + grid(i,j).location(1), bottom*ones(N,1)];
        w = w_h;
    elseif face == 3
        r = [right*ones(N,1), flip(y)*dy/2 + grid(i,j).location(2)];
        w = w_v;
    elseif face == 4
        r = [flip(x)*dx/2 + grid(i,j).location(1), top*ones(N,1)];
        w = w_h;
    end
    
elseif strcmp(quad_method,'Trapz') == 1
    if face == 1
        r = [left*ones(N,1), transpose(linspace(bottom,top,N))];
    elseif face == 2
        r = [transpose(linspace(left,right,N)), bottom*ones(N,1)];
    elseif face == 3
        r = [right*ones(N,1), transpose(linspace(bottom,top,N))];
    elseif face == 4
        r = [transpose(linspace(left,right,N)), top*ones(N,1)];
    end
end

