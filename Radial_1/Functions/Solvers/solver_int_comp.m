function zn = solver_int_comp(grid,i,j,source_gp,r_s,N,N_sub,quad_rule)

    zn = grid(i,j);
    
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
    
    x_bounds = [zn.location(1) - zn.edge_lengths(1)/2, zn.location(1) + zn.edge_lengths(1)/2];
    y_bounds = [zn.location(2) - zn.edge_lengths(2)/2, zn.location(2) + zn.edge_lengths(2)/2];
    [zn.avg_flux,zn.avg_N_absorb] = interior_comp(grid,i,j,r_s,source_gp,nodes,weights,x_bounds,y_bounds,N_sub,trapz_flag);
    
    %If the source lies in the center of a zone, slightly offset it to 
    %prevent opt_depth from equaling 0 for plotting
    if isequal(r_s,zn.location)
        R_cent = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2) + 1/20*grid(i,j).edge_lengths(1); 
    else
        R_cent = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2);
    end
    zn.sigma_bar_zone = average_xs(zn.location(1:2),r_s,source_gp,grid,i,j);
    zn.opt_depth_zone = R_cent*zn.sigma_bar_zone;
        
    
    function [avg_flux,avg_absorb] = interior_comp(grid,i,j,r_s,source_gp,nodes,weights,x_bounds,y_bounds,N_sub,trapz_flag)
        dim = grid(end,end).location + grid(end,end).edge_lengths/2;
        x_ints = linspace(x_bounds(1),x_bounds(2),N_sub+1);
        y_ints = linspace(y_bounds(1),y_bounds(2),N_sub+1);
        sub_flux = zeros(N_sub);
        sub_absorb = zeros(N_sub);

        for m = 1:length(x_ints)-1
            dx = x_ints(m+1) - x_ints(m);
            sub_nodes_x = 1/2*(1+nodes)*dx + x_ints(m);
            sub_weights_x = 1/2*weights*dx;
            for n = 1:length(y_ints)-1
                dy = y_ints(n+1) - y_ints(n);
                sub_nodes_y = 1/2*(1+nodes)*dy + y_ints(n);
                sub_weights_y = 1/2*weights*dy;
                [X,Y] = meshgrid(sub_nodes_x,sub_nodes_y);
                [WX,WY] = meshgrid(sub_weights_x,sub_weights_y);
                flux_eval = zeros(length(sub_nodes_x));
                absorb_eval = zeros(length(sub_nodes_x));
                
                for k = 1:length(sub_nodes_x)
                    for p = 1:length(sub_nodes_y)
                        r = [X(k,p),Y(k,p)];
                        R = norm(r-r_s);
                        sigma_bar = average_xs(r,r_s,source_gp,grid,i,j);
                        flux_eval(k,p) = exp(-sigma_bar*R)/R;
                        absorb_eval(k,p) = flux_eval(k,p)*xs_quad_points(grid(i,j),r,dim)*dx*dy;
                    end
                end
                
                if trapz_flag == 0
                    sub_flux(m,n) = sum(sum(WX.*WY.*flux_eval))/(dx*dy);
                    sub_absorb(m,n) = sum(sum(WX.*WY.*absorb_eval))/(dx*dy);
                else
                    sub_flux(m,n) = trapz(Y(:,1),trapz(X(1,:),flux_eval,2))/(dx*dy);
                    sub_absorb(m,n) = trapz(Y(:,1),trapz(X(1,:),absorb_eval,2))/(dx*dy);
                end
            end
        end
            
        avg_flux = sum(sum(sub_flux*dx*dy))/((x_bounds(2)-x_bounds(1))*(y_bounds(2)-y_bounds(1)));
        avg_absorb = sum(sum(sub_absorb));
    end
end

