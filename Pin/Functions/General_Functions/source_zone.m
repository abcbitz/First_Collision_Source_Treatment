function gp = source_zone(r_s,grid)
%gp is the gridpoint of the zone containing the source
%Ff the source lies on the grid, gp will contain an array of all zones
%adjacent to the source - 2 if on a line, 4 if on a corner. If along an
%edge, it will return zones outside of grid (they will not be used later
%on).

num_zones = size(grid);
gp = [];

for i = 1:num_zones(1)
    for j = 1:num_zones(2)
        left = grid(i,j).location(1) - grid(i,j).edge_lengths(1)/2;
        right = grid(i,j).location(1) + grid(i,j).edge_lengths(1)/2;
        bottom = grid(i,j).location(2) - grid(i,j).edge_lengths(2)/2;
        top = grid(i,j).location(2) + grid(i,j).edge_lengths(2)/2;
        if r_s(1) >= left && r_s(1) <= right && r_s(2) >= bottom && r_s(2) <= top
            gp(1,:) = [i,j];
            if r_s(1) == left 
                gp(size(gp,1)+1,:) = [i-1,j];
            elseif r_s(1) == right
                gp(size(gp,1)+1,:) = [i+1,j];
            end
            
            if r_s(2) == bottom 
                gp(size(gp,1)+1,:) = [i,j-1];
            elseif r_s(2) == top
                gp(size(gp,1)+1,:) = [i,j+1];
            end
            
            if size(gp,1) == 3
                %if in a corner along the edges, it may not be in top right
                if r_s(1) == 0.0 && r_s(2) == 0.0
                    gp(size(gp,1)+1,:) = [i-1,j-1];
                elseif r_s(1) == 0.0
                    gp(size(gp,1)+1,:) = [i-1,j+1];
                elseif r_s(2) == 0.0
                    gp(size(gp,1)+1,:) = [i+1,j-1];
                else
                    gp(size(gp,1)+1,:) = [i+1,j+1];
                end
            end
            
            return
        end
    end
end


