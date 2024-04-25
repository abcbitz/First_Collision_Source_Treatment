function zn = solver_surf_comp(grid,i,j,source_gp,r_s,N,N_sub,quad_rule)

    zn = grid(i,j);
    
    dx = zn.edge_lengths(1);
    dy = zn.edge_lengths(2);
    
    trapz_flag = 0;
    if strcmp(quad_rule,'Gauss')
        [nodes,weights] = lgwt(N,1,-1);
        nodes = transpose(nodes);
        weights = transpose(-weights);
        
    elseif strcmp(quad_rule,'Gauss-Lobatto')
        [nodes,weights] = lglnodes(N-1);
        nodes = transpose(flip(nodes));
        weights = transpose(weights);
        
    else
        nodes = linspace(-1,1,N);
        weights = ones(1,N);
        trapz_flag = 1;
    end
    
    %(x,y) is the top right corner of the zone (the far edges)
    %boudns is zone bounds in x or y
    
    %Face 1
    if i == 1
        x = zn.location(1) - zn.edge_lengths(1)/2;
        bounds = [zn.location(2) - zn.edge_lengths(2)/2, zn.location(2) + zn.edge_lengths(2)/2];
        zn.total_current(1) = gauss_kron(grid,i,j,x,[-1,0],r_s,source_gp,nodes,weights,bounds,N_sub,trapz_flag);
    else
        zn.total_current(1) = -grid(i-1,j).total_current(3);
    end
    
    %Face 2
    if j == 1
        y = zn.location(2) - zn.edge_lengths(2)/2;
        bounds = [zn.location(1) - zn.edge_lengths(1)/2, zn.location(1) + zn.edge_lengths(1)/2];
        zn.total_current(2) = gauss_kron(grid,i,j,y,[0,-1],r_s,source_gp,nodes,weights,bounds,N_sub,trapz_flag);
    else 
        zn.total_current(2) = -grid(i,j-1).total_current(4);
    end
    
    %Face 3
    x = zn.location(1) + zn.edge_lengths(1)/2;
    bounds = [zn.location(2) - zn.edge_lengths(2)/2, zn.location(2) + zn.edge_lengths(2)/2];
    zn.total_current(3) = gauss_kron(grid,i,j,x,[1,0],r_s,source_gp,nodes,weights,bounds,N_sub,trapz_flag);
    
    %Face 4
    y = zn.location(2) + zn.edge_lengths(2)/2;
    bounds = [zn.location(1) - zn.edge_lengths(1)/2, zn.location(1) + zn.edge_lengths(1)/2];
    zn.total_current(4) = gauss_kron(grid,i,j,y,[0,1],r_s,source_gp,nodes,weights,bounds,N_sub,trapz_flag);
    
    if ismember([i,j],source_gp,'rows')
        zn.avg_N_absorb = 2*pi/size(source_gp,1) - sum(zn.total_current.*[dy,dx,dy,dx]);
    else
        zn.avg_N_absorb = sum(-zn.total_current.*[dy,dx,dy,dx]);
    end
    
    global use_func;
    if use_func == 1
        dim = grid(end,end).location + grid(end,end).edge_lengths/2;
        xs_avg = avg_zonal_xs(zn,dim);
        zn.avg_flux = zn.avg_N_absorb/(prod(zn.edge_lengths)*xs_avg);
    else
        zn.avg_flux = zn.avg_N_absorb/(prod(zn.edge_lengths)*zn.sigma_a);
    end
    
    %If the source lies in the center of a zone, slightly offset it to 
    %prevent opt_depth from equaling 0 for plotting
    if isequal(r_s,zn.location)
        R_cent = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2) + 1/20*dx; 
    else
        R_cent = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2);
    end
    zn.sigma_bar_zone = average_xs(zn.location(1:2),r_s,source_gp,grid,i,j);
    zn.opt_depth_zone = R_cent*zn.sigma_bar_zone;
        
    
    function total_current = gauss_kron(grid,i,j,xy,n,r_s,source_gp,nodes,weights,bounds,N_sub,trapz_flag)
        ints = transpose(linspace(bounds(1),bounds(2),N_sub+1));
        approx = zeros(N_sub,1);

        for m = 1:length(ints)-1
            sub_nodes = 1/2*(1+nodes)*(ints(m+1) - ints(m)) + ints(m);
            sub_weights = 1/2*weights*(ints(m+1) - ints(m));
            f_eval = zeros(1,length(sub_nodes));
            for k = 1:length(sub_nodes)
                if n(2) == 0
                    r = [xy,sub_nodes(k)];
                else 
                    r = [sub_nodes(k),xy];
                end
                R = sqrt(sum((r-r_s).^2));
                sigma_bar = average_xs(r,r_s,source_gp,grid,i,j);
                f_eval(k) = exp(-sigma_bar*R)/R * dot(n,(r-r_s)/R);
                
                if i == 2 && j == 3 && k == 1 && n(2) == 1
                    fprintf('\nX_bounds = %.4f %.4f\n',ints(m),ints(m+1))
                end
                if i == 2 && j == 3 && n(2) == 1
                    fprintf('%.4f  %.8f\n',sub_nodes(k),f_eval(k))
                end
            end
            
            if trapz_flag == 0
                approx(m) = sum(sub_weights.*f_eval)/(ints(m+1) - ints(m));
            else
                approx(m) = trapz(sub_nodes,f_eval)/(ints(m+1) - ints(m));
            end
            
            if i == 2 && j == 3 && n(2) == 1
                fprintf('sub_current = %.12f\n',approx(m))
            end
        end

        total_current = sum(approx.*(ints(2:end)-ints(1:end-1)))/(bounds(2)-bounds(1));
        
        if i == 2 && j == 3 && n(2) == 1
            fprintf('\nTotal_current = %.12f\n\n',total_current)
        end
    end
end

