function xs_avg = avg_zonal_xs(zn,dim)
%Finds the average xs across a zone. Used for the surface method when
%integrating the xs function instead of performing the ray trace on the
%mesh.

N = 5;

[nodes,weights] = lgwt(N,1,-1);
nodes = transpose(nodes);
weights = transpose(-weights);
x_bounds = [zn.location(1) - zn.edge_lengths(1)/2, zn.location(1) + zn.edge_lengths(1)/2];
y_bounds = [zn.location(2) - zn.edge_lengths(2)/2, zn.location(2) + zn.edge_lengths(2)/2];

dx = zn.edge_lengths(1);
nodes_x = 1/2*(1+nodes)*dx + x_bounds(1);
weights_x = 1/2*weights*dx;

dy = zn.edge_lengths(2);
nodes_y = 1/2*(1+nodes)*dy + y_bounds(1);
weights_y = 1/2*weights*dy;

[X,Y] = meshgrid(nodes_x,nodes_y);
[WX,WY] = meshgrid(weights_x,weights_y);

xs = zeros(N);
for i = 1:N
    for j = 1:N
        r = [X(i,j),Y(i,j)];
        xs(i,j) = xs_generator(r,dim);
    end
end

xs_avg = sum(sum(WX.*WY.*xs))/(dx*dy);
        
        
        