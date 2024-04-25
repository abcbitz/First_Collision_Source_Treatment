function zn = solver_surf(grid,i,j,source_gp,r_s,N,quad_rule,neg_fix)

    zn = grid(i,j);

    dx = grid(i,j).edge_lengths(1);
    dy = grid(i,j).edge_lengths(2);

    %for debugging
    %if i == 8 && j == 1 
    %    fprintf('')
    %end
    
    %would run faster if intead of calling quad_points_surf everytime you
    %just adjusted the points from -1 to 1 to fit the zone

    %Left face (1)
    if i == 1
        [r,w] = quad_points_surf(grid,i,j,1,N,quad_rule);
        for n = 1:N
            R = sqrt(sum((r(n,:) - r_s).^2));
            zn.sigma_bar_1(n,:) = [r(n,1), r(n,2), average_xs(r(n,:),r_s,source_gp,grid,i,j)];
            zn.opt_depth_1(n,:) = [r(n,1), r(n,2), zn.sigma_bar_1(n,3)*R];
            zn.currents_1(n,:) = [r(n,1), r(n,2), exp(-zn.opt_depth_1(n,3))/R * dot([-1,0],(r(n,:)-r_s)/R)];
        end

        if strcmp(quad_rule,'Gauss') || strcmp(quad_rule,'Gauss-Lobatto') == 1
            zn.total_current(1) = sum(w.*zn.currents_1(:,3))/dy; %Gauss quadrature
        elseif strcmp(quad_rule,'Trapz') == 1
            zn.total_current(1) = trapz(r(:,2),zn.currents_1(:,3))/dy; %Trapezoidal rule
        end

        if ismember(r_s,r,'rows') && size(source_gp,1) == 4 %excludes case where on a single gridline
            zn.total_current(1) = 0.0;
        end

    else
        zn.sigma_bar_1 = grid(i-1,j).sigma_bar_3;
        zn.opt_depth_1 = grid(i-1,j).opt_depth_3;
        zn.currents_1(:,1:2) = grid(i-1,j).currents_3(:,1:2);
        zn.currents_1(:,3) = -grid(i-1,j).currents_3(:,3);
        zn.total_current(1) = -grid(i-1,j).total_current(3);
    end


    %Bottom face (2)
    if j == 1
        [r,w] = quad_points_surf(grid,i,j,2,N,quad_rule);
        for n = 1:N
            R = sqrt(sum((r(n,:) - r_s).^2));
            zn.sigma_bar_2(n,:) = [r(n,1), r(n,2), average_xs(r(n,:),r_s,source_gp,grid,i,j)];
            zn.opt_depth_2(n,:) = [r(n,1), r(n,2), zn.sigma_bar_2(n,3)*R];
            zn.currents_2(n,:) = [r(n,1), r(n,2), exp(-zn.opt_depth_2(n,3))/R * dot([0,-1],(r(n,:)-r_s)/R)];
        end

        if strcmp(quad_rule,'Gauss') || strcmp(quad_rule,'Gauss-Lobatto') == 1
            zn.total_current(2) = sum(w.*zn.currents_2(:,3))/dx; %Gauss quadrature
        elseif strcmp(quad_rule,'Trapz') == 1
            zn.total_current(2) = trapz(r(:,1),zn.currents_2(:,3))/dx; %Trapezoidal rule
        end

        if ismember(r_s,r,'rows') && size(source_gp,1) == 4 %excludes case where on a single gridline
            zn.total_current(2) = 0.0;
        end

    else
        zn.sigma_bar_2 = grid(i,j-1).sigma_bar_4;
        zn.opt_depth_2 = grid(i,j-1).opt_depth_4;
        zn.currents_2(:,1:2) = grid(i,j-1).currents_4(:,1:2);
        zn.currents_2(:,3) = -grid(i,j-1).currents_4(:,3);
        zn.total_current(2) = -grid(i,j-1).total_current(4);
    end


    %Right face (3)
    [r,w] = quad_points_surf(grid,i,j,3,N,quad_rule);
    for n = 1:N
        R = sqrt(sum((r(n,:) - r_s).^2));
        zn.sigma_bar_3(n,:) = [r(n,1), r(n,2), average_xs(r(n,:),r_s,source_gp,grid,i,j)];
        zn.opt_depth_3(n,:) = [r(n,1), r(n,2), zn.sigma_bar_3(n,3)*R];
        zn.currents_3(n,:) = [r(n,1), r(n,2), exp(-zn.opt_depth_3(n,3))/R * dot([1,0],(r(n,:)-r_s)/R)];
    end

    if strcmp(quad_rule,'Gauss') || strcmp(quad_rule,'Gauss-Lobatto') == 1
        zn.total_current(3) = sum(w.*zn.currents_3(:,3))/dy; %Gauss quadrature
    elseif strcmp(quad_rule,'Trapz') == 1
        zn.total_current(3) = trapz(r(:,2),zn.currents_3(:,3))/dy; %Trapezoidal rule
    end

    if ismember(r_s,r,'rows') && size(source_gp,1) == 4 %excludes case where on a single gridline
        zn.total_current(3) = 0.0;
    end


    %Top face (4)
    [r,w] = quad_points_surf(grid,i,j,4,N,quad_rule);
    for n = 1:N
        R = sqrt(sum((r(n,:) - r_s).^2));
        zn.sigma_bar_4(n,:) = [r(n,1), r(n,2), average_xs(r(n,:),r_s,source_gp,grid,i,j)];
        zn.opt_depth_4(n,:) = [r(n,1), r(n,2), zn.sigma_bar_4(n,3)*R];
        zn.currents_4(n,:) = [r(n,1), r(n,2), exp(-zn.opt_depth_4(n,3))/R * dot([0,1],(r(n,:)-r_s)/R)];
    end

    if strcmp(quad_rule,'Gauss') || strcmp(quad_rule,'Gauss-Lobatto') == 1
        zn.total_current(4) = sum(w.*zn.currents_4(:,3))/dx; %Gauss quadrature
    elseif strcmp(quad_rule,'Trapz') == 1
        zn.total_current(4) = trapz(r(:,1),zn.currents_4(:,3))/dx; %Trapezoidal rule
    end

    if ismember(r_s,r,'rows') && size(source_gp,1) == 4 %excludes case where on a single gridline
        zn.total_current(4) = 0.0;
    end


    %If the source lies in the center of a zone, slightly offset it to 
    %prevent opt_depth from equaling 0 for plotting    
    if isequal(r_s,zn.location)
        R = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2) + 1/20*dx;
    else
        R = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2);
    end
    zn.sigma_bar_zone = average_xs(zn.location(1:2),r_s,source_gp,grid,i,j);
    zn.opt_depth_zone = R*zn.sigma_bar_zone;
    
    %Calculate number of absorptions
    if ismember([i,j],source_gp,'rows')
        zn.avg_N_absorb = 2*pi/size(source_gp,1) - sum(zn.total_current.*[dy,dx,dy,dx]);
    else
        zn.avg_N_absorb = sum(-zn.total_current.*[dy,dx,dy,dx]); %abs/s
    end
    
    global use_func;
    %If using the function, need to take an average across the zone
    %Should zn.sigma_a just be the zonal avg instead of the midpoint?
    if use_func == 1
        dim = grid(end,end).location + grid(end,end).edge_lengths/2;
        xs_avg = avg_zonal_xs(zn,dim);
        zn.avg_flux = zn.avg_N_absorb/(prod(zn.edge_lengths)*xs_avg);
    else
        zn.avg_flux = zn.avg_N_absorb/(prod(zn.edge_lengths)*zn.sigma_a);
    end

    %In comp, might want to have a flag telling it to recalculate the
    %current across faces 1 and 2 if being called by neg_fix
    if zn.avg_N_absorb < 0.0
        fprintf('Negative flux at (%d,%d)\n',i,j);  %might get erased
        if neg_fix == 1
            zn = solver_surf_comp(grid,i,j,source_gp,r_s,N,6);
        end
    end
end


