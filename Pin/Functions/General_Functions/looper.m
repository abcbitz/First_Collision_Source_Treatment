function grid = looper(num_zones,dim,N,q,r_s,solver,quad_rule,neg_fix,N_sub,max_error)
%would like to add the solver specific variables to two cell arrays and
%read them in that way

grid = build_grid(num_zones,dim);
s = fprintf('Percent Done: 0.0%%\n');

source_gp = source_zone(r_s,grid);

for i = 1:num_zones(1)
    for j = 1:num_zones(2)
        %call solvers
        if strcmp(solver,'solver_surf')
            grid(i,j) = solver_surf(grid,i,j,source_gp,r_s,N,quad_rule,neg_fix);
        elseif strcmp(solver,'solver_surf_comp')
            grid(i,j) = solver_surf_comp(grid,i,j,source_gp,r_s,N,N_sub,quad_rule);
        elseif strcmp(solver,'solver_surf_adapt')
            grid(i,j) = solver_surf_adapt(grid,i,j,source_gp,r_s,N,N_sub,max_error);
        elseif strcmp(solver,'solver_int')
            grid(i,j) = solver_int(grid,i,j,source_gp,r_s,N,quad_rule);
        elseif strcmp(solver,'solver_int_comp')
            grid(i,j) = solver_int_comp(grid,i,j,source_gp,r_s,N,N_sub,quad_rule);
        end
        
    end
    
    if rem(i,5) == 0.0
        s = fprintf('Percent Done: %%%.2f\n',i/num_zones(1)*100);
    end
%     if i == num_zones(1)
%         fprintf(repmat('\b',1,s))
%     end
end

%Scale everything by q at end to preserve precision
for i = 1:num_zones(1)
    for j = 1:num_zones(2)
        grid(i,j).flux = grid(i,j).flux*q;
        grid(i,j).N_absorb = grid(i,j).N_absorb*q;
        grid(i,j).sub_currents_1 = grid(i,j).sub_currents_1*q;
        grid(i,j).sub_currents_2 = grid(i,j).sub_currents_2*q;
        grid(i,j).sub_currents_3 = grid(i,j).sub_currents_3*q;
        grid(i,j).sub_currents_4 = grid(i,j).sub_currents_4*q;
        if strcmp(solver,'solver_surf')
            grid(i,j).currents_1(:,3) = grid(i,j).currents_1(:,3)*q;
            grid(i,j).currents_2(:,3) = grid(i,j).currents_2(:,3)*q;
            grid(i,j).currents_3(:,3) = grid(i,j).currents_3(:,3)*q;
            grid(i,j).currents_4(:,3) = grid(i,j).currents_4(:,3)*q;
        end
        grid(i,j).total_current = grid(i,j).total_current*q;
        grid(i,j).avg_N_absorb = grid(i,j).avg_N_absorb*q;
        grid(i,j).avg_flux = grid(i,j).avg_flux*q;
    end
end

disp('')

