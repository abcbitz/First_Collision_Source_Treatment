function extra_axes(num_zones,flag)

return

%Turn off figure box - remove all ticks
box off;
ax1 = gca;
%Add box back - with lower x-axis ticks
ax2 = axes('position',ax1.Position,'box','on','ytick',[],'xtick',[],'color','none');

%Add second x and y axes
%1D slice
if flag == 1
    ax3 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','xlim',[0,num_zones(1)],'YLim',ax1.YLim,'YScale',ax1.YScale);
%Temperature plot
else
    ax3 = axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','xlim',[0,num_zones(1)],'YLim',[0,num_zones(2)]);
end

