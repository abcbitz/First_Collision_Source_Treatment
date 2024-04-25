function [absorb_error_tot,absorb_error_avg,absorb_error_med,absorb_error_L1,...
    absorb_error_L2,absorb_error_L1_rem,absorb_error_L2_rem,leakage_error,...
    cons_error,wp_rel_error,wp_L1] = ...
    plotter_error_dist(grid,grid_true,q,r_s,dim,N,solver,includePlots,savePlots,exclude,wp)

%input vars
set_num = 5;
quad_meth = 'Gauss';
%quad_meth = 'Gauss-Lobatto';
extra_plots = 0;
lineout_err = 0;
rows = [2,8];
columns = [2,8];
bounds = [1e-16,1e-8];


edge_lengths = [grid(1,1).edge_lengths(2),grid(1,1).edge_lengths(1),grid(1,1).edge_lengths(2),grid(1,1).edge_lengths(1)];
source_gp = source_zone(r_s,grid);
num_zones = size(grid);

x = zeros(num_zones);
y = zeros(num_zones);
X = zeros(num_zones);
Y = zeros(num_zones);
flux = zeros(num_zones);
flux_true = zeros(num_zones);
J_in = zeros(num_zones);
J_out = zeros(num_zones);
N_in = zeros(num_zones);
N_out = zeros(num_zones);
N_abs = zeros(num_zones);
N_abs_true = zeros(num_zones);
xs = zeros(num_zones);
vol = zeros(num_zones);
opt_depth = zeros(num_zones);
sigma_bar = zeros(num_zones);
tot_leakage = 0.0;
tot_leakage_true = 0.0;
leakage_L1 = 0.0;
currents = zeros(num_zones(1)+1,num_zones(2)+1,4);
currents_true = zeros(num_zones(1)+1,num_zones(2)+1,4);

for i = 1:num_zones(1)
    for j = 1:num_zones(2)
        x(i,j) = grid(i,j).location(1);
        y(i,j) = grid(i,j).location(2);
        X(i,j) = grid(i,j).location(1) - grid(i,j).edge_lengths(1)/2;
        Y(i,j) = grid(i,j).location(2) - grid(i,j).edge_lengths(2)/2;
        vol(i,j) = prod(grid(i,j).edge_lengths);
        flux(i,j) = grid(i,j).avg_flux;
        flux_true(i,j) = grid_true(i,j).avg_flux;
        N_abs(i,j) = grid(i,j).avg_N_absorb;
        N_abs_true(i,j) = grid_true(i,j).avg_N_absorb;
        xs(i,j) = grid(i,j).sigma_a;
        opt_depth(i,j) = grid(i,j).opt_depth_zone;
        sigma_bar(i,j) = grid(i,j).sigma_bar_zone; 
        
        if contains(solver,'surf')
            currents(i,j,:) = grid(i,j).total_current;
        else
            currents(i,j,:) = grid_true(i,j).total_current;
        end
        currents_true(i,j,:) = grid_true(i,j).total_current;
        
        %Compute current in and current out
        for k = 1:4
            if currents(i,j,k) < 0.0
                J_in(i,j) = J_in(i,j) + abs(currents(i,j,k)); 
                N_in(i,j) = N_in(i,j) + abs(currents(i,j,k)*edge_lengths(k));
            else
                J_out(i,j) = J_out(i,j) + abs(currents(i,j,k)); 
                N_out(i,j) = N_out(i,j) + abs(currents(i,j,k)*edge_lengths(k));
            end
        end
        
        %Compute current in for source zones
        if ismember([i,j],source_gp,'rows')
            J_in(i,j) = 1;  %This is just a placeholder since current doesn't apply here
            N_in(i,j) = 2*pi/size(source_gp,1)*q;
        end
        
        %Calculate leakage leaving the problem   
        %Current is always positive when it is leaking from a zone
        % - check current on outer edges
        if i == 1
            tot_leakage = tot_leakage + currents(i,j,1)*grid(i,j).edge_lengths(2);
            tot_leakage_true = tot_leakage_true + grid_true(i,j).total_current(1)*grid_true(i,j).edge_lengths(2);
            leakage_L1 = leakage_L1 + abs(currents(i,j,1)-grid_true(i,j).total_current(1))*grid(i,j).edge_lengths(2);
        end
        if j == 1
            tot_leakage = tot_leakage + currents(i,j,2)*grid(i,j).edge_lengths(1);
            tot_leakage_true = tot_leakage_true + grid_true(i,j).total_current(2)*grid_true(i,j).edge_lengths(1);
            leakage_L1 = leakage_L1 + abs(currents(i,j,2)-grid_true(i,j).total_current(2))*grid(i,j).edge_lengths(1);
        end
        if i == num_zones(1)
            tot_leakage = tot_leakage + currents(i,j,3)*grid(i,j).edge_lengths(2);
            tot_leakage_true = tot_leakage_true + grid_true(i,j).total_current(3)*grid_true(i,j).edge_lengths(2);
            leakage_L1 = leakage_L1 + abs(currents(i,j,3)-grid_true(i,j).total_current(3))*grid(i,j).edge_lengths(2);
        end
        if j == num_zones(2)
            tot_leakage = tot_leakage + currents(i,j,4)*grid(i,j).edge_lengths(1);
            tot_leakage_true = tot_leakage_true + grid_true(i,j).total_current(4)*grid_true(i,j).edge_lengths(1);
            leakage_L1 = leakage_L1 + abs(currents(i,j,4)-grid_true(i,j).total_current(4))*grid(i,j).edge_lengths(1);
        end
    end
end
leakage_L1 = leakage_L1/tot_leakage_true;


if wp == -1
    wp_rel_error = 0;
    wp_L1 = 0;
else
    dx = dim(1)/wp;
    x_wp = zeros(wp,1);
    for i = 1:wp
        x_wp(i) = (i-0.5)*dx;
    end
    dy = dim(2)/wp;
    y_wp = zeros(wp,1);
    for i = 1:wp
        y_wp(i) = (i-0.5)*dy;
    end
    wp_rel_error = zeros(wp);
    N_abs_wp = zeros(wp);
    N_abs_true_wp = zeros(wp);
    for p = 1:wp
        for k = 1:wp
            for i = 1:num_zones(1)
                for j = 1:num_zones(2)
                    if abs(x(i,j)-x_wp(p)) <= 1e-10 && abs(y(i,j)-y_wp(k)) <= 1e-10
                        N_abs_wp(p,k) = N_abs(i,j);
                        N_abs_true_wp(p,k) = N_abs_true(i,j);
                        wp_rel_error(p,k) = abs(N_abs(i,j)-N_abs_true(i,j))/N_abs_true(i,j);
                    end
                end
            end
        end
    end
    wp_L1 = sum(sum(abs(N_abs_wp-N_abs_true_wp)))/sum(sum(N_abs_true_wp));
end

if isequal(r_s,[0.0,0.0]) || isequal(r_s,[dim(1),0.0]) || isequal(r_s,[0.0,dim(2)]) || isequal(r_s,dim) %corners
    S = pi/2*q;
elseif r_s(1) == 0.0 || r_s(1) == dim(1) || r_s(2) == 0.0 || r_s(2) == dim(2) %edges
    S = pi*q;
else
    S = pi*2*q; %inside grid
end

total_abs = sum(sum(N_abs));
fprintf('Calculated Total Absorption: %0.6e\n', total_abs);
true_total_abs = sum(sum(N_abs_true));
fprintf('True Total Absorption: %0.6e\n', true_total_abs)
absorb_error_tot = abs(true_total_abs-total_abs)/true_total_abs;
fprintf('Total Absorption Error: %0.6e\n', absorb_error_tot)
absorb_error_avg = mean(mean(abs(N_abs-N_abs_true)./N_abs_true));
fprintf('Mean Absorption Error: %0.6e\n', absorb_error_avg)
absorb_error_med = median(median(abs(N_abs-N_abs_true)./N_abs_true));
fprintf('Median Absorption Error: %0.6e\n', absorb_error_med)
absorb_error_L1 = sum(sum(abs(N_abs-N_abs_true)))/sum(sum(N_abs_true));
fprintf('L1 Relative Absorb Error: %0.6e\n', absorb_error_L1)
absorb_error_L2 = sqrt(sum(sum((N_abs-N_abs_true).^2)))/sum(sum(N_abs_true));  
fprintf('L2 Relative Absorb Error: %0.6e\n\n', absorb_error_L2)

if ~exist('exclude','var')  %percentage of zones to exclude in bottom left corner
    exclude = 10;
end

top = zeros(2,1);
bot = 0.0;
for i = ceil(exclude/100*num_zones(1)):num_zones(1)
    for j = ceil(exclude/100*num_zones(2)):num_zones(2)
        top(1) = top(1) + abs(N_abs(i,j)-N_abs_true(i,j));
        top(2) = top(2) + (N_abs(i,j)-N_abs_true(i,j))^2;
        bot = bot + N_abs_true(i,j);
    end
end
absorb_error_L1_rem = top(1)/bot;
absorb_error_L2_rem = sqrt(top(2))/bot;

fprintf('Removing zones (i,j) = [1:%d,1:%d]\n',ceil(num_zones(1)*0.2),ceil(num_zones(2)*0.2))
fprintf('L1 Relative Absorb Error: %0.6e\n', absorb_error_L1_rem)
fprintf('L2 Relative Absorb Error: %0.6e\n\n', absorb_error_L2_rem)

%Leakge error for the interior point method will equal 0
%tot_leakage_true = S - true_total_abs;
fprintf('Calculated Leakage: %0.15e\n', tot_leakage);
fprintf('True Leakage: %0.15e\n', tot_leakage_true)
leakage_error = abs(tot_leakage_true-tot_leakage)/tot_leakage_true;
fprintf('Total Leakage Error: %0.6e\n\n', leakage_error)
fprintf('L1 Relative Leakage Error: %0.6e\n\n', leakage_L1)

%Print out conservation values
fprintf('Conservation: \nSource Rate (S): %0.8e\n', S);
fprintf('Total Absorption Rate (N_abs): %0.8e\n', total_abs);
fprintf('Leakage (J_edge): %0.8e\n', tot_leakage);
fprintf('N_src - N_abs - N_leakage = %0.6e\n', S - total_abs - tot_leakage);
cons_error = abs(S - total_abs - tot_leakage)/S;
fprintf('(N_src - N_abs - N_leakage)/N_src = %0.6e\n', cons_error);
fprintf('Conservation Acheived - ');
if abs(S- total_abs - tot_leakage)/S <= 1.0e-6
    fprintf('yes\n\n')
else
    fprintf('no\n\n')
end


figsize = [50 50 800 600];
fontsize = 15;
gridlines = 1;
n = 200;

%leaving this in so the move in the main doesnt error
folder = sprintf('.\\Runs\\%s\\XY_%i_%i__N_%i',solver,num_zones(1),num_zones(2),N);
mkdir(folder)
save(sprintf('%s\\data',folder))


folder = sprintf('.\\Runs\\Set_%i\\error_dist\\%s\\%s',set_num,quad_meth,solver);
mkdir(folder)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HEAT MAPS
%Things get messed with a bit for visual clarity. 

%Other than X and Y, the actual value here doesnt matter, just place
%holders to get the last zone to show up because MATLAB sucks
for i = 1:num_zones(1)
    X(i,num_zones(2)+1) = grid(i,num_zones(2)).location(1) - grid(i,num_zones(2)).edge_lengths(1)/2;
    Y(i,num_zones(2)+1) = grid(i,num_zones(2)).location(2) + grid(i,num_zones(2)).edge_lengths(2)/2;
    x(i,num_zones(2)+1) = grid(i,num_zones(2)).location(1);
    y(i,num_zones(2)+1) = grid(i,num_zones(2)).location(2) + grid(i,num_zones(2)).edge_lengths(2);
    opt_depth(i,num_zones(2)+1) = opt_depth(i,num_zones(2));
    flux(i,num_zones(2)+1) = flux(i,num_zones(2));
    flux_true(i,num_zones(2)+1) = flux_true(i,num_zones(2));
    N_abs(i,num_zones(2)+1) = N_abs(i,num_zones(2));
    N_abs_true(i,num_zones(2)+1) = N_abs_true(i,num_zones(2));
    J_in(i,num_zones(2)+1) = J_in(i,num_zones(2));
    J_out(i,num_zones(2)+1) = J_out(i,num_zones(2));    
    N_in(i,num_zones(2)+1) = N_in(i,num_zones(2));
    N_out(i,num_zones(2)+1) = N_out(i,num_zones(2));
    xs(i,num_zones(2)+1) = xs(i,num_zones(2));
    sigma_bar(i,num_zones(2)+1) = sigma_bar(i,num_zones(2));
    vol(i,num_zones(2)+1) = vol(i,num_zones(2));
end
for j = 1:num_zones(2)
    X(num_zones(1)+1,j) = grid(num_zones(1),j).location(1) + grid(num_zones(1),j).edge_lengths(1)/2;
    Y(num_zones(1)+1,j) = grid(num_zones(1),j).location(2) - grid(num_zones(1),j).edge_lengths(2)/2;
    x(num_zones(1)+1,j) = grid(num_zones(1),j).location(1) + grid(num_zones(1),j).edge_lengths(1);
    y(num_zones(1)+1,j) = grid(num_zones(1),j).location(2);
    opt_depth(num_zones(1)+1,j) = opt_depth(num_zones(1),j);
    flux(num_zones(1)+1,j) = flux(num_zones(1),j);
    flux_true(num_zones(1)+1,j) = flux_true(num_zones(1),j);
    N_abs(num_zones(1)+1,j) = N_abs(num_zones(1),j);
    N_abs_true(num_zones(1)+1,j) = N_abs_true(num_zones(1),j);
    J_in(num_zones(1)+1,j) = J_in(num_zones(1),j);
    J_out(num_zones(1)+1,j) = J_out(num_zones(1),j);
    N_in(num_zones(1)+1,j) = N_in(num_zones(1),j);
    N_out(num_zones(1)+1,j) = N_out(num_zones(1),j);
    xs(num_zones(1)+1,j) = xs(num_zones(1),j);
    sigma_bar(num_zones(1)+1,j) = sigma_bar(num_zones(1),j);
    vol(num_zones(1)+1,j) = vol(num_zones(1),j);
end
X(num_zones(1)+1,num_zones(2)+1) = X(num_zones(1)+1,num_zones(2));
Y(num_zones(1)+1,num_zones(2)+1) = Y(num_zones(1),num_zones(2)+1);
x(num_zones(1)+1,num_zones(2)+1) = x(num_zones(1)+1,num_zones(2));
y(num_zones(1)+1,num_zones(2)+1) = y(num_zones(1),num_zones(2)+1);
opt_depth(num_zones(1)+1,num_zones(2)+1) = opt_depth(num_zones(1)+1,num_zones(2));
flux(num_zones(1)+1,num_zones(2)+1) = flux(num_zones(1)+1,num_zones(2));
flux_true(num_zones(1)+1,num_zones(2)+1) = flux_true(num_zones(1)+1,num_zones(2));
N_abs(num_zones(1)+1,num_zones(2)+1) = N_abs(num_zones(1)+1,num_zones(2));
N_abs_true(num_zones(1)+1,num_zones(2)+1) = N_abs_true(num_zones(1)+1,num_zones(2));
J_in(num_zones(1)+1,num_zones(2)+1) = J_in(num_zones(1)+1,num_zones(2));
J_out(num_zones(1)+1,num_zones(2)+1) = J_out(num_zones(1)+1,num_zones(2));
N_in(num_zones(1)+1,num_zones(2)+1) = J_in(num_zones(1)+1,num_zones(2));
N_out(num_zones(1)+1,num_zones(2)+1) = J_out(num_zones(1)+1,num_zones(2));
xs(num_zones(1)+1,num_zones(2)+1) = xs(num_zones(1)+1,num_zones(2));
sigma_bar(num_zones(1)+1,num_zones(2)+1) = sigma_bar(num_zones(1)+1,num_zones(2));
vol(num_zones(1)+1,num_zones(2)+1) = vol(num_zones(1)+1,num_zones(2));


figure(n)
clf(n)
set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
rel_error = abs(N_abs - N_abs_true)./N_abs_true;
s = pcolor(X,Y,rel_error);
if gridlines ~= 1
    s.EdgeColor = 'none';
end
title('Relative Absorption Error',' ','Fontsize',18)
xlabel('x [cm]')
ylabel('y [cm]')
colormap turbo
colorbar('Location','westoutside')
set(gca,'ColorScale','log')
mm = [10^floor(log10(min(min(rel_error)))),10^ceil(log10(max(max(rel_error))))];
%caxis(mm)
caxis([1e-14,1e-1])
extra_axes(num_zones,2)
%save_plot(figure(n),sprintf('%s\\rel_abs_error',folder),savePlots)
save_plot(figure(n),sprintf('%s\\rel_abs_error__%ix%i__%i_rays',folder,num_zones(1),num_zones(2),N),1)
n = n + 1;


current_error = abs(currents-currents_true)./abs(currents_true);
if extra_plots == 1
    figure(n)
    clf(n)
    set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
    rel_error_R = rel_error.*sqrt(x.^2+y.^2);
    s = pcolor(X,Y,rel_error_R);
    if gridlines ~= 1
        s.EdgeColor = 'none';
    end
    title('Relative Absorption Error x R',' ','Fontsize',18)
    xlabel('x [cm]')
    ylabel('y [cm]')
    colormap turbo
    colorbar('Location','westoutside')
    set(gca,'ColorScale','log')
    mm = [10^floor(log10(min(min(rel_error_R)))),10^ceil(log10(max(max(rel_error_R))))];
    %caxis(mm)
    caxis([1e-14,1e-1])
    extra_axes(num_zones,2)
    save_plot(figure(n),sprintf('%s\\rel_abs_error_R__%ix%i__%i_rays',folder,num_zones(1),num_zones(2),N),1)
    n = n + 1;
    
    figure(n)
    clf(n)
    set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
    abs_error = abs(N_abs - N_abs_true);
    s = pcolor(X,Y,abs_error);
    if gridlines ~= 1
        s.EdgeColor = 'none';
    end
    title('Absorption Error',' ','Fontsize',18)
    xlabel('x [cm]')
    ylabel('y [cm]')
    colormap turbo
    colorbar('Location','westoutside')
    set(gca,'ColorScale','log')
    mm = [10^floor(log10(min(min(abs_error)))),10^ceil(log10(max(max(abs_error))))];
    %caxis(mm)
    caxis([1e-14,1e-1])
    extra_axes(num_zones,2)
    save_plot(figure(n),sprintf('%s\\abs_error__%ix%i__%i_rays',folder,num_zones(1),num_zones(2),N),1)
    %save_plot(figure(n),sprintf('.\\Runs\\Set_13\\Run__Gauss\\abs_error_XY_%i__%i_rays__surf',num_zones(1),N),1)
    n = n + 1;
    
    figure(n)
    clf(n)
    set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
    rel_flux_error = abs(flux - flux_true)./flux_true;
    s = pcolor(X,Y,rel_flux_error);
    if gridlines ~= 1
        s.EdgeColor = 'none';
    end
    title('Relative Flux Error',' ','Fontsize',18)
    xlabel('x [cm]')
    ylabel('y [cm]')
    colormap turbo
    colorbar('Location','westoutside')
    set(gca,'ColorScale','log')
    mm = [10^floor(log10(min(min(rel_flux_error)))),10^ceil(log10(max(max(rel_flux_error))))];
    %caxis(mm)
    caxis([1e-14,1e-1])
    extra_axes(num_zones,2)
    save_plot(figure(n),sprintf('%s\\rel_flux_error__%ix%i__%i_rays',folder,num_zones(1),num_zones(2),N),1)
    n = n + 1;
    
    for i = 1:4
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        s = pcolor(X,Y,current_error(:,:,i));
        if gridlines ~= 1
            s.EdgeColor = 'none';
        end
        title(sprintf('Current Error %i',i),' ','Fontsize',18)
        xlabel('x [cm]')
        ylabel('y [cm]')
        colormap turbo
        colorbar('Location','westoutside')
        set(gca,'ColorScale','log')
        caxis([1e-14,1e-1])
        extra_axes(num_zones,2)
        save_plot(figure(n),sprintf('%s\\current_error_%i__%ix%i__%i_rays',folder,i,num_zones(1),num_zones(2),N),1)
        n = n+1;
    end
end

if lineout_err == 1
    for j = rows
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        semilogy(X(1:end-1,j),rel_error(1:end-1,j),'x-')
        title(sprintf('Relative Absorbtion Error - j = %i', j),' ','Fontsize',18)
        xlabel('x [cm]')
        ylabel('Rel Absorp Error')
        ylim(bounds)
        extra_axes(num_zones,1)
        save_plot(figure(n),sprintf('%s\\rel_abs_error_j%i__%ix%i__%i_rays',folder,j,num_zones(1),num_zones(2),N),1)
        n = n+1;
        
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        semilogy(X(1:end-1,j),current_error(1:end-1,j,3),'x-')
        title(sprintf('Relative Current 3 Error - j = %i', j),' ','Fontsize',18)
        xlabel('x [cm]')
        ylabel('Rel Cur Error')
        ylim(bounds)
        extra_axes(num_zones,1)
        save_plot(figure(n),sprintf('%s\\rel_cur3_error_j%i__%ix%i__%i_rays',folder,j,num_zones(1),num_zones(2),N),1)
        n = n+1;
        
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        semilogy(X(1:end-1,j),current_error(1:end-1,j,4),'x-')
        title(sprintf('Relative Current 4 Error - j = %i', j),' ','Fontsize',18)
        xlabel('x [cm]')
        ylabel('Rel Cur Error')
        ylim(bounds)
        extra_axes(num_zones,1)
        save_plot(figure(n),sprintf('%s\\rel_cur4_error_j%i__%ix%i__%i_rays',folder,j,num_zones(1),num_zones(2),N),1)
        n = n+1;
    end
    
    for i = columns
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        semilogy(Y(i,1:end-1),rel_error(i,1:end-1),'x-')
        title(sprintf('Relative Absorbtion Error - i = %i', i),' ','Fontsize',18)
        xlabel('y [cm]')
        ylabel('Rel Absorp Error')
        ylim(bounds)
        extra_axes(num_zones,1)
        save_plot(figure(n),sprintf('%s\\rel_abs_error_i%i__%ix%i__%i_rays',folder,i,num_zones(1),num_zones(2),N),1)
        n = n+1;
        
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        semilogy(Y(i,1:end-1),current_error(i,1:end-1,3),'x-')
        title(sprintf('Relative Current 3 Error - i = %i', i),' ','Fontsize',18)
        xlabel('y [cm]')
        ylabel('Rel Cur Error')
        ylim(bounds)
        extra_axes(num_zones,1)
        save_plot(figure(n),sprintf('%s\\rel_cur3_error_i%i__%ix%i__%i_rays',folder,i,num_zones(1),num_zones(2),N),1)
        n = n+1;
        
        figure(n)
        clf(n)
        set(gcf,'Position', figsize,'DefaultAxesFontSize',fontsize)
        semilogy(Y(i,1:end-1),current_error(i,1:end-1,4),'x-')
        title(sprintf('Relative Current 4 Error - i = %i', i),' ','Fontsize',18)
        xlabel('y [cm]')
        ylabel('Rel Cur Error')
        ylim(bounds)
        extra_axes(num_zones,1)
        save_plot(figure(n),sprintf('%s\\rel_cur4_error_i%i__%ix%i__%i_rays',folder,i,num_zones(1),num_zones(2),N),1)
        n = n+1;
    end
end


