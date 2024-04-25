classdef zone
    
    properties
        gridpoint       %location in grid (1,1),(1,2),(2,1)...
        location        %zone center location
        edge_lengths    %width, height
        sigma_a         %absorption cross section
        
        %interior point
        WX
        WY
        X
        Y
        xs              %xs at quad points
        sigma_bar       %[angle, sigma_bar]
        opt_depth       %[angle, optical depth] or optical depth
        flux
        N_absorb
        
        %surface method
        sigma_bar_1     %[x, y, sigma_bar]
        sigma_bar_2     %[x, y, sigma_bar]
        sigma_bar_3     %[x, y, sigma_bar]
        sigma_bar_4     %[x, y, sigma_bar]
        opt_depth_1     %[x, y, opt_depth]
        opt_depth_2     %[x, y, opt_depth]
        opt_depth_3     %[x, y, opt_depth]
        opt_depth_4     %[x, y, opt_depth]
        currents_1      %[x, y, angular flux]
        currents_2      %[x, y, angular flux]
        currents_3      %[x, y, angular flux]
        currents_4      %[x, y, angular flux]
        total_current         %current
        
        %Both - zonal values
        sigma_bar_zone
        opt_depth_zone
        avg_N_absorb        %number of particles absorbed
        avg_flux        %uncollided flux
        
        %true solution
        intervals_1
        intervals_2
        intervals_3
        intervals_4
        sub_currents_1
        sub_currents_2
        sub_currents_3
        sub_currents_4
        num_divs
    end
    
    
    methods
        function obj = zone(gridpoint,edge_lengths,location,sigma_a)
            obj.gridpoint = gridpoint;
            obj.edge_lengths = edge_lengths;
            obj.location = location;
            obj.sigma_a = sigma_a; 
        end
    end
end
