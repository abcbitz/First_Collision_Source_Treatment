function zn = solver_surf_adapt(grid,i,j,source_gp,r_s,N,max_divs,max_error)

    if N == 3
        nodes = [-0.7745966692414833770359,0.0,0.7745966692414833770359];
        weights = [0.5555555555555555555556,0.8888888888888888888889,0.5555555555555555555556];

    elseif N == 5
        nodes = [-0.9258200997725514615666,-0.5773502691896257645092,0.0,...
            0.5773502691896257645092,0.9258200997725514615666];
        weights = [0.197979797979797979798,0.4909090909090909090909,0.6222222222222222222222,...
            0.4909090909090909090909,0.197979797979797979798];

    elseif N == 7
        nodes = [-0.9604912687080202834235,-0.7745966692414833770359,-0.4342437493468025580021,...
            0.0,0.4342437493468025580021,0.7745966692414833770359,0.9604912687080202834235];
        weights = [0.1046562260264672651938,0.2684880898683334407286,...
            0.4013974147759622229051,0.4509165386584741423451,0.4013974147759622229051,...
            0.2684880898683334407286,0.1046562260264672651938];

    elseif N == 9
        nodes = [-0.9765602507375731115345,-0.861136311594052575224,-0.6402862174963099824047,...
            -0.3399810435848562648027,0,0.3399810435848562648027,0.6402862174963099824047,...
            0.861136311594052575224,0.9765602507375731115345];
        weights = [0.06297737366547301476549,0.1700536053357227268027,0.2667983404522844480328,...
            0.326949189601451629558,0.346442981890136361681,0.326949189601451629558,0.2667983404522844480328,...
            0.1700536053357227268027,0.06297737366547301476549];

    elseif N == 11
        nodes = [-0.9840853600948424644962,-0.9061798459386639927976,-0.7541667265708492204408,...
            -0.5384693101056830910363,-0.2796304131617831934135,0,0.2796304131617831934135,...
            0.5384693101056830910363,0.7541667265708492204408,0.9061798459386639927976,...
            0.9840853600948424644962];
        weights = [0.04258203675108183286451,0.1152333166224733940246,0.186800796556492657468,...
            0.2410403392286475866999,0.272849801912558922341,0.2829874178574912132043,...
            0.272849801912558922341,0.2410403392286475866999,0.186800796556492657468,...
            0.1152333166224733940246,0.04258203675108183286451];

    elseif N == 13
        nodes = [-0.9887032026126788575047,-0.9324695142031520278123,-0.8213733408650279400457,...
            -0.6612093864662645136614,-0.4631182124753046121568,-0.2386191860831969086305,0,...
            0.2386191860831969086305,0.4631182124753046121568,0.6612093864662645136614,...
            0.8213733408650279400457,0.9324695142031520278123,0.9887032026126788575047];
        weights = [0.030396154119819768852,0.0836944404469066261328,0.1373206046344469230872,...
            0.181071994323137615187,0.213209652271962279163,0.233770864116994406623,...
            0.2410725801734647619106,0.233770864116994406623,0.213209652271962279163,...
            0.181071994323137615187,0.1373206046344469230872,0.0836944404469066261328,...
            0.030396154119819768852];

    elseif N == 15
        nodes = [-0.991455371120813,-0.949107912342759,-0.864864423359769,...
            -0.741531185599394,-0.586087235467691,-0.405845151377397,...
            -0.207784955007898,0.000000000000000,0.207784955007898,...
            0.405845151377397,0.586087235467691,0.741531185599394,...
            0.864864423359769,0.949107912342759,0.991455371120813];
        weights = [0.022935322010529,0.063092092629979,0.104790010322250,...
            0.140653259715525,0.169004726639267,0.190350578064785,...
            0.204432940075298,0.209482141084728,0.204432940075298,...
            0.190350578064785,0.169004726639267,0.140653259715525,...
            0.104790010322250,0.063092092629979,0.022935322010529];
    end
    
    zn = grid(i,j);
    
    dx = zn.edge_lengths(1);
    dy = zn.edge_lengths(2);
    
    %(x,y) is the top right corner of the zone (the far edges)
    %bounds is zone bounds in x or y
    
    %Face 1
    if i == 1
        x = zn.location(1) - zn.edge_lengths(1)/2;
        bounds = [zn.location(2) - zn.edge_lengths(2)/2, zn.location(2) + zn.edge_lengths(2)/2];
        %[zn.total_current(1),zn.intervals_1,zn.sub_currents_1,zn.num_divs(1)] = gauss_kron(grid,i,j,x,[-1,0],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
        [zn.total_current(1),~,~,zn.num_divs(1)] = gauss_kron(grid,i,j,x,[-1,0],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
    else
        zn.num_divs(1) = grid(i-1,j).num_divs(3);
        zn.total_current(1) = -grid(i-1,j).total_current(3);
    end
    
    %Face 2
    if j == 1
        y = zn.location(2) - zn.edge_lengths(2)/2;
        bounds = [zn.location(1) - zn.edge_lengths(1)/2, zn.location(1) + zn.edge_lengths(1)/2];
        %[zn.total_current(2),zn.intervals_2,zn.sub_currents_2,zn.num_divs(2)] = gauss_kron(grid,i,j,y,[0,-1],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
        [zn.total_current(2),~,~,zn.num_divs(2)] = gauss_kron(grid,i,j,y,[0,-1],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
    else 
        zn.num_divs(2) = grid(i,j-1).num_divs(4);
        zn.total_current(2) = -grid(i,j-1).total_current(4);
    end
    
    %Face 3
    x = zn.location(1) + zn.edge_lengths(1)/2;
    bounds = [zn.location(2) - zn.edge_lengths(2)/2, zn.location(2) + zn.edge_lengths(2)/2];
    [zn.total_current(3),zn.intervals_3,zn.sub_currents_3,zn.num_divs(3)] = gauss_kron(grid,i,j,x,[1,0],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
    %[zn.total_current(3),~,~,zn.num_divs(3)] = gauss_kron(grid,i,j,x,[1,0],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
    
    %Face 4
    y = zn.location(2) + zn.edge_lengths(2)/2;
    bounds = [zn.location(1) - zn.edge_lengths(1)/2, zn.location(1) + zn.edge_lengths(1)/2];
    [zn.total_current(4),zn.intervals_4,zn.sub_currents_4,zn.num_divs(4)] = gauss_kron(grid,i,j,y,[0,1],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
    %[zn.total_current(4),~,~,zn.num_divs(4)] = gauss_kron(grid,i,j,y,[0,1],r_s,source_gp,nodes,weights,0,bounds,0.0,[],0.0,0,max_error,max_divs);
    
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

        
    
    function [total_current,ints,approx,divs] = gauss_kron(grid,i,j,xy,n,r_s,source_gp,nodes,...
            weights,bad_ints,ints,approx,prev_ints,prev_approx,divs,max_er,max_divs)
        %ints is updated to contain each of the different sub-intervals
        %approx is the solution for those intervals
        %bad_ints is the list of intervals that still need refining

        for m = 1:length(ints)-1
            if bad_ints(m) == 0
                sub_nodes = 1/2*(1+nodes)*(ints(m+1) - ints(m)) + ints(m);
                sub_weights = 1/2*weights*(ints(m+1) - ints(m));
                f_eval = zeros(1,length(sub_nodes));
                for k = 1:length(sub_nodes)
                    if n(2) == 0 %if n is horiz, d is the x coord
                        r = [xy,sub_nodes(k)];
                    else %if n is vert, d is the y coord
                        r = [sub_nodes(k),xy];
                    end
                    R = sqrt(sum((r-r_s).^2));
                    sigma_bar = average_xs(r,r_s,source_gp,grid,i,j);
                    f_eval(k) = exp(-sigma_bar*R)/R * dot(n,(r-r_s)/R);
                end
                approx(m) = sum(sub_weights.*f_eval);
            end
        end
        
        total_current = sum(approx)/(ints(end)-ints(1));
        prev_total_current = sum(prev_approx)/(ints(end)-ints(1));
        
        %if the current intervals solution has converged, add it to approx and continue to next interval
        er = abs((total_current - prev_total_current)/(total_current+1e-17));   %add 1e-17 to avoid dividing by 0
        if i == 1 && j == 4 && n(2) == 1
            disp('')
        end
        if length(ints) > 2 && (er <= max_er || divs >= max_divs)
%             if er >= max_er && divs >= max_divs
%                 fprintf('Function not converging correctly (%.5e) at zone %d,%d - n = [%d,%d]\n',er,i,j,n(1),n(2))
%                 fprintf('')
%             end
            return
        else
            new_ints = ints(1);
            new_approx = approx(1);
            bad_ints = bad_ints*0;  %technically unnecessary
            if length(ints) == 2    %initial interval
                new_ints = [ints(1), ints(1) + 1/2*(ints(2)-ints(1)), ints(2)];
                new_approx = [0,0];
                bad_ints = [0,0];
            else
                count1 = 2; %ints counter
                count2 = 2; %new_ints counter
                %Identify bad sub-intervals and sub-divide them
                for p = 2:length(prev_ints)
                    %no sub-division previously, so already know it converged
                    if prev_ints(p) == ints(count1)
                        new_ints(count2) = ints(count1);
                        new_approx(count2-1) = approx(count1-1);
                        bad_ints(count2-1) = 1;
                        count1 = count1 + 1;
                        count2 = count2 + 1;
                        
                    %sub-division occured last iteration
                    else
                        er = abs(prev_approx(p-1)-approx(count1-1)-approx(count1))/abs(approx(count1-1)+approx(count1)+1e-17);
                        if er >= max_er
                            new_ints(count2:count2+3) = [ints(count1)-1/2*(ints(count1)-ints(count1-1)), ints(count1), ints(count1)+1/2*(ints(count1+1)-ints(count1)), ints(count1+1)];
                            new_approx(count2-1:count2+2) = [0,0,0,0];
                            bad_ints(count2-1:count2+2) = [0,0,0,0];
                            count1 = count1 + 2;
                            count2 = count2 + 4;
                        else
                            new_ints(count2:count2+1) = [ints(count1),ints(count1+1)];
                            new_approx(count2-1:count2) = [approx(count1-1),approx(count1)];
                            bad_ints(count2-1:count2) = [1,1];
                            count1 = count1 + 2;
                            count2 = count2 + 2;
                        end
                    end
                end
            end
            
            divs = divs + 1;
            [total_current,ints,approx,divs] = gauss_kron(grid,i,j,xy,n,r_s,source_gp,nodes,weights,bad_ints,new_ints,new_approx,ints,approx,divs,max_er,max_divs);
        end
    end
end

