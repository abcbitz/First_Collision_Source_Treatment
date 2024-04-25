function grid = build_grid(num_zones,dim)

grid = repmat(zone([1,1],1,1,1),num_zones(1),num_zones(2));
edge_lengths = dim./num_zones;

for i=1:num_zones(1)
    for j=1:num_zones(2)
        r = [edge_lengths(1)*(i-1/2),edge_lengths(2)*(j-1/2)];
        xs = xs_generator(r,dim);
        %xs = xs_generator(flip(r),dim);
        grid(i,j) = zone([i,j], edge_lengths, r, xs);
    end
end



