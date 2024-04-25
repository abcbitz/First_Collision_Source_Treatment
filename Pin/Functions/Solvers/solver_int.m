function zn = solver_int(grid,i,j,source_gp,r_s,N,quad_rule)

    zn = grid(i,j);

    dx = grid(i,j).edge_lengths(1);
    dy = grid(i,j).edge_lengths(2);
    
    dim = grid(end,end).location + grid(end,end).edge_lengths/2;

    %for pausing in debug
    %if i == 8 && j == 1 
    %    fprintf('')
    %end

    %Check and see if a quad point lies on the source. If so, try
    %shifting the quadrature by increasing order by 1. If this does not
    %work, the source is on a gridpoint, so switch to Gauss. 

    [zn.X, zn.Y, zn.WX, zn.WY] = quad_points_cent(grid,i,j,N,quad_rule);

    %check if any point is on the source
    if ismember(r_s,reshape(cat(2,zn.X',zn.Y'),[],2),'rows')
        N = N + 1;
        [zn.X, zn.Y, zn.WX, zn.WY] = quad_points_cent(grid,i,j,N,quad_rule);

        %if still on source, switch to Gauss
        if ismember(r_s,reshape(cat(2,zn.X',zn.Y'),[],2),'rows')
            quad_rule = 'Gauss';
            N = N - 1;
            [zn.X, zn.Y, zn.WX, zn.WY] = quad_points_cent(grid,i,j,N,quad_rule);
        end
    end
    
    %Loop over the quad points and compute the scalar flux at each
    for m = 1:N
        for n = 1:N
            r = [zn.X(m,n), zn.Y(m,n)];
            R = sqrt(sum((r-r_s).^2));
            zn.sigma_bar(m,n) = average_xs(r,r_s,source_gp,grid,i,j);
            zn.opt_depth(m,n) = zn.sigma_bar(m,n)*R;
            zn.flux(m,n) = exp(-zn.opt_depth(m,n))/R;
            zn.xs(m,n) = xs_quad_points(zn,r,dim);
            zn.N_absorb(m,n) = zn.flux(m,n)*zn.xs(m,n)*dx*dy; 
        end
    end

    %integrate across the quad points to get the avg flux
    if strcmp(quad_rule,'Gauss') == 1 || strcmp(quad_rule,'Gauss-Lobatto') == 1 %Gauss quadrature
        zn.avg_flux = sum(sum(zn.WX.*zn.WY.*zn.flux))/(dx*dy);
        zn.avg_N_absorb = sum(sum(zn.WX.*zn.WY.*zn.N_absorb))/(dx*dy);   
    elseif strcmp(quad_rule,'Trapz') == 1 %Trapezoidal rule
        zn.avg_flux = trapz(zn.Y(:,1),trapz(zn.X(1,:),zn.flux,2))/(dx*dy); 
        zn.avg_N_absorb = trapz(zn.Y(:,1),trapz(zn.X(1,:),zn.N_absorb,2))/(dx*dy);  
    end

    %Zone center values for plotting
    if isequal(source_gp,[i,j]) && isequal(r_s,zn.location)
        R = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2) + 1/20*dx; %To stop opt_depth from equaling 0 for plotting
        zn.sigma_bar_zone = zn.sigma_a;
        zn.opt_depth_zone = R*zn.sigma_bar_zone;
    else
        R = sqrt((zn.location(1)-r_s(1))^2 + (zn.location(2)-r_s(2))^2);
        zn.sigma_bar_zone = average_xs(zn.location(1:2),r_s,source_gp,grid,i,j);
        zn.opt_depth_zone = R*zn.sigma_bar_zone;
    end
end


