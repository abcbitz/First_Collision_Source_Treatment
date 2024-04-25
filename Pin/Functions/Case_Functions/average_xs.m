function sigma_bar = average_xs(r,r_s,source_gp,grid,i_end,j_end)

    %fprintf('(i, j, x, y) = [%i, %i, %.4f, %.4f] \n',i_end,j_end,r(1),r(2))
%     if i_end == 2 && j_end == 1
%         fprintf('(i, j, x, y) = [%i, %i, %.4f, %.4f] \n',i_end,j_end,r(1),r(2))
%         disp('')
%     end

    global use_func;
    if use_func == 1
        %dont remember how this works but it does...
        g = @(t) (5* ((1-t)*r_s(1) + t*r(1)).^2 + 1) * sqrt((r(1)-r_s(1))^2 + (r(2)-r_s(2))^2);
        sigma_bar = integral(g,0,1,'RelTol',1e-15)/norm(r-r_s);
        return
    end


    %upper right one, unless it is outside the grid
    %as long as i,j_src are in the grid, im not sure if it matters which zone
    %is used...
    x_vals = unique(source_gp(:,1));
    if x_vals(end) <= size(grid,1)
        i_src = x_vals(end);
    else
        i_src = x_vals(1);
    end

    y_vals = unique(source_gp(:,2));
    if y_vals(end) <= size(grid,2)
        j_src = y_vals(end);
    else
        j_src = y_vals(1);
    end

    R_tot = norm(r - r_s);
    
    x_prev = r_s(1);
    y_prev = r_s(2);
    i = i_src;
    j = j_src;
    R = zeros(abs(i_end - source_gp(1,1)) + abs(j_end - source_gp(1,2)) + 1, 5);
    sigma = zeros(abs(i_end - source_gp(1,1)) + abs(j_end - source_gp(1,2)) + 1, 1);

    if r(1) > r_s(1) && r(2) > r_s(2) %from 0 to 90 deg
        x_c = grid(i_src,j_src).location(1) + grid(i_src,j_src).edge_lengths(1)/2;  %x coord corner
        y_c = grid(i_src,j_src).location(2) + grid(i_src,j_src).edge_lengths(2)/2;  %y coord corner
        x = x_func(r,r_s,y_c);  %where the ray crosses the next horizontal face
        y = y_func(r,r_s,x_c);  %where the ray crosses the next vertical face
        count = 1;
        while true
            sigma(count,1) = grid(i,j).sigma_a;
            if i == i_end && j == j_end
                R(count,:) = [x_prev,y_prev,r(1),r(2),sqrt((r(1)-x_prev)^2 + (r(2)-y_prev)^2)];
                break
                
            elseif abs(y-y_c) <= 1e-15    %passes through corner
                R(count,:) = [x_prev,y_prev,x_c,y_c,sqrt((x_c-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x_c;
                y_prev = y_c;
                i = i + 1;
                j = j + 1;
                if i > i_end || j > j_end   %if r is on a corner, needs special conditional
                    break
                end
                x_c = grid(i,j).location(1) + grid(i,j).edge_lengths(1)/2;
                y_c = grid(i,j).location(2) + grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
                y = y_func(r,r_s,x_c);

            elseif y < y_c  %passes through vertical face
                R(count,:) = [x_prev,y_prev,x_c,y,sqrt((x_c-x_prev)^2 + (y-y_prev)^2)];
                x_prev = x_c;   %where the ray is currently in x
                y_prev = y;     %where the ray is currently in y
                i = i + 1;
                x_c = grid(i,j).location(1) + grid(i,j).edge_lengths(1)/2;  %the x coord of the top right corner
                y = y_func(r,r_s,x_c);  %the value of y for the ray at x_c

            elseif y > y_c   %passes through horizontal face
                R(count,:) = [x_prev,y_prev,x,y_c,sqrt((x-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x;
                y_prev = y_c;
                j = j + 1;
                y_c = grid(i,j).location(2) + grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
            end

            count = count + 1;
        end


    elseif r(1) < r_s(1) && r(2) > r_s(2)   %90 to 180 deg
        x_c = grid(i_src,j_src).location(1) - grid(i_src,j_src).edge_lengths(1)/2;
        y_c = grid(i_src,j_src).location(2) + grid(i_src,j_src).edge_lengths(2)/2;
        x = x_func(r,r_s,y_c);
        y = y_func(r,r_s,x_c);
        count = 1;
        while true
            sigma(count,1) = grid(i,j).sigma_a;
            if i == i_end && j == j_end
                R(count,:) = [x_prev,y_prev,r(1),r(2),sqrt((r(1)-x_prev)^2 + (r(2)-y_prev)^2)];
                break
                
            elseif abs(y-y_c) <= 1e-15
                R(count,:) = [x_prev,y_prev,x_c,y_c,sqrt((x_c-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x_c;
                y_prev = y_c;
                i = i - 1;
                j = j + 1;
                if i < i_end || j > j_end   %if r is on a corner, needs special conditional
                    break
                end
                x_c = grid(i,j).location(1) - grid(i,j).edge_lengths(1)/2;
                y_c = grid(i,j).location(2) + grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
                y = y_func(r,r_s,x_c);

            elseif y < y_c
                R(count,:) = [x_prev,y_prev,x_c,y,sqrt((x_c-x_prev)^2 + (y-y_prev)^2)];
                x_prev = x_c;
                y_prev = y;
                i = i - 1;
                x_c = grid(i,j).location(1) - grid(i,j).edge_lengths(1)/2;
                y = y_func(r,r_s,x_c);

            elseif y > y_c
                R(count,:) = [x_prev,y_prev,x,y_c,sqrt((x-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x;
                y_prev = y_c;
                j = j + 1;
                y_c = grid(i,j).location(2) + grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
            end

            count = count + 1;
        end


    elseif r(1) < r_s(1) && r(2) < r_s(2)   %180 to 270 deg
        x_c = grid(i_src,j_src).location(1) - grid(i_src,j_src).edge_lengths(1)/2;
        y_c = grid(i_src,j_src).location(2) - grid(i_src,j_src).edge_lengths(2)/2;
        x = x_func(r,r_s,y_c);
        y = y_func(r,r_s,x_c);
        count = 1;
        while true
            sigma(count,1) = grid(i,j).sigma_a;
            if i == i_end && j == j_end
                R(count,:) = [x_prev,y_prev,r(1),r(2),sqrt((r(1)-x_prev)^2 + (r(2)-y_prev)^2)];
                break
                
            elseif abs(y-y_c) <= 1e-15
                R(count,:) = [x_prev,y_prev,x_c,y_c,sqrt((x_c-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x_c;
                y_prev = y_c;
                i = i - 1;
                j = j - 1;
                if i < i_end || j < j_end   %if r is on a corner, needs special conditional
                    break
                end
                x_c = grid(i,j).location(1) - grid(i,j).edge_lengths(1)/2;
                y_c = grid(i,j).location(2) - grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
                y = y_func(r,r_s,x_c);

            elseif y > y_c
                R(count,:) = [x_prev,y_prev,x_c,y,sqrt((x_c-x_prev)^2 + (y-y_prev)^2)];            
                x_prev = x_c;
                y_prev = y;
                i = i - 1;
                x_c = grid(i,j).location(1) - grid(i,j).edge_lengths(1)/2;
                y = y_func(r,r_s,x_c);

            elseif y < y_c
                R(count,:) = [x_prev,y_prev,x,y_c,sqrt((x-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x;
                y_prev = y_c;
                j = j - 1;
                y_c = grid(i,j).location(2) - grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
            end

            count = count + 1;
        end


    elseif r(1) > r_s(1) && r(2) < r_s(2)   %270 to 360 deg
        x_c = grid(i_src,j_src).location(1) + grid(i_src,j_src).edge_lengths(1)/2;
        y_c = grid(i_src,j_src).location(2) - grid(i_src,j_src).edge_lengths(2)/2;
        x = x_func(r,r_s,y_c);
        y = y_func(r,r_s,x_c);
        count = 1;
        while true
            sigma(count,1) = grid(i,j).sigma_a;
            if i == i_end && j == j_end
                R(count,:) = [x_prev,y_prev,r(1),r(2),sqrt((r(1)-x_prev)^2 + (r(2)-y_prev)^2)];
                break

            elseif abs(y-y_c) <= 1e-15
                R(count,:) = [x_prev,y_prev,x_c,y_c,sqrt((x_c-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x_c;
                y_prev = y_c;
                i = i + 1;
                j = j - 1;
                if i > i_end || j < j_end   %if r is on a corner, needs special conditional
                    break
                end
                x_c = grid(i,j).location(1) + grid(i,j).edge_lengths(1)/2;
                y_c = grid(i,j).location(2) - grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
                y = y_func(r,r_s,x_c);
                
            elseif y > y_c
                R(count,:) = [x_prev,y_prev,x_c,y,sqrt((x_c-x_prev)^2 + (y-y_prev)^2)]; 
                x_prev = x_c;
                y_prev = y;
                i = i + 1;
                x_c = grid(i,j).location(1) + grid(i,j).edge_lengths(1)/2;
                y = y_func(r,r_s,x_c);

            elseif y < y_c
                R(count,:) = [x_prev,y_prev,x,y_c,sqrt((x-x_prev)^2 + (y_c-y_prev)^2)];
                x_prev = x;
                y_prev = y_c;
                j = j - 1;
                y_c = grid(i,j).location(2) - grid(i,j).edge_lengths(2)/2;
                x = x_func(r,r_s,y_c);
            end
            
            count = count + 1;
        end


    elseif r(2) == r_s(2)   %0 or 180 degrees
        %What xs do you use when the source is in line with a face?
        %Really need to make sure the right source indices are being used
        %What do these look like when at a node? 
        % - if going right, should be fine. If left, I think you need the src
        % to be on the right face, not the left. - R(1) will equal 0, but this
        % might be fine...
        count = 2;
        if r(1) > r_s(1)    %going right
            R(1,5) = grid(i_src,j_src).location(1) + grid(i_src,j_src).edge_lengths(1)/2 - r_s(1);    %src to right edge
            sigma(1,1) = grid(i_src,j_src).sigma_a;
            for i = (i_src+1):(i_end-1)
                R(count,5) = grid(i,j_src).edge_lengths(1);
                sigma(count,1) = grid(i,j_src).sigma_a;
                count = count + 1;
            end
            R(count,5) = r(1) - (grid(i_end,j_end).location(1) - grid(i_end,j_end).edge_lengths(1)/2);
            sigma(count,1) = grid(i_end,j_end).sigma_a;

        else    %going left
            R(1,5) = r_s(1) - grid(i_src,j_src).location(1) + grid(i_src,j_src).edge_lengths(1)/2;    %src to left edge
            sigma(1,1) = grid(i_src,j_src).sigma_a;
            for i = (i_src-1):-1:(i_end+1)
                R(count,5) = grid(i,j_src).edge_lengths(1);
                sigma(count,1) = grid(i,j_src).sigma_a;
                count = count + 1;
            end
            R(count,5) = grid(i_end,j_end).location(1) + grid(i_end,j_end).edge_lengths(1)/2 - r(1);
            sigma(count,1) = grid(i_end,j_end).sigma_a;
        end


    elseif r(1) == r_s(1)   %90 or 270 degrees
        count = 2;
        if r(2) > r_s(2)    %going up
            R(1,5) = grid(i_src,j_src).location(2) + grid(i_src,j_src).edge_lengths(2)/2 - r_s(2); 
            sigma(1,1) = grid(i_src,j_src).sigma_a;
            for j = (j_src+1):(j_end-1)
                R(count,5) = grid(i_src,j).edge_lengths(2);
                sigma(count,1) = grid(i_src,j).sigma_a;
                count = count + 1;
            end
            R(count,5) = r(2) - (grid(i_end,j_end).location(2) - grid(i_end,j_end).edge_lengths(2)/2);
            sigma(count,1) = grid(i_end,j_end).sigma_a;

        else    %going down
            R(1,5) = r_s(2)  - grid(i_src,j_src).location(2) + grid(i_src,j_src).edge_lengths(2)/2;
            sigma(1,1) = grid(i_src,j_src).sigma_a;
            for j = (j_src-1):-1:(j_end+1)
                R(count,5) = grid(i_src,j).edge_lengths(2);
                sigma(count,1) = grid(i_src,j).sigma_a;
                count = count + 1;
            end
            R(count,5) = grid(i_end,j_end).location(2) + grid(i_end,j_end).edge_lengths(2)/2 - r(2);
            sigma(count,1) = grid(i_end,j_end).sigma_a;
        end

    else
        fprintf('Something is wrong at %i,%i \n',i_end,j_end)
    end

    sigma_bar = sum(sigma.*R(:,5))/sum(R(:,5));


    %Rounding these functions because precision error is causing issues
    %in the conditionals. Shouldn't have an effect.
    function x = x_func(r,r_s,y)
        x = (r(1)-r_s(1))/(r(2)-r_s(2))*(y-r_s(2)) + r_s(1);
        %x = round((r(1)-r_s(1))/(r(2)-r_s(2))*(y-r_s(2)) + r_s(1),14);
    end

    function y = y_func(r,r_s,x)
        y = (r(2)-r_s(2))/(r(1)-r_s(1))*(x-r_s(1)) + r_s(2);
        %y = round((r(2)-r_s(2))/(r(1)-r_s(1))*(x-r_s(1)) + r_s(2),14);
    end
end
