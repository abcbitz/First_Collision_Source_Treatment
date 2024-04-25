function save_plot(fig,name,savePlots)

if savePlots == 1
    saveas(fig,name,'png')
    saveas(fig,name,'fig')
end

