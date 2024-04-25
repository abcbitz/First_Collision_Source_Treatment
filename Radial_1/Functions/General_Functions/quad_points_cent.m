function [X,Y,WX,WY] = quad_points_cent(grid,i,j,N,quad_method)

dx = grid(i,j).edge_lengths(1);
dy = grid(i,j).edge_lengths(2);

left = grid(i,j).location(1) - dx/2;
right = grid(i,j).location(1) + dx/2;
bottom = grid(i,j).location(2) - dy/2;
top = grid(i,j).location(2) + dy/2;

w = [];
r = [];

WX = [];
WY = [];
X = [];
Y = [];

if strcmp(quad_method,'Gauss') == 1
    [x,wx] = lgwt(N,left,right);
    [y,wy] = lgwt(N,bottom,top);
    [WX,WY] = meshgrid(wx,wy);
    w = reshape(cat(2,WX',WY'),[],2);
    [X,Y] = meshgrid(x,y);
    r = reshape(cat(2,X',Y'),[],2);
    
elseif strcmp(quad_method,'Gauss-Lobatto') == 1
    [p,w,~] = lglnodes(N-1);
    [WX,WY] = meshgrid(w*dx/2,w*dy/2);
    w = reshape(cat(2,WX',WY'),[],2);
    x = p*dx/2 + grid(i,j).location(1);
    y = p*dy/2 + grid(i,j).location(2);
    [X,Y] = meshgrid(x,y);
    r = reshape(cat(2,X',Y'),[],2);
    
elseif strcmp(quad_method,'Trapz') == 1
    y = transpose(linspace(bottom,top,N));
    x = transpose(linspace(left,right,N));
    [X,Y] = meshgrid(x,y);
    r = reshape(cat(2,X',Y'),[],2);
end

