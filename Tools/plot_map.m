function plot_map(lonT, latT, quantity, C_label, Fontsize, Fontsize_Ticks)
% function that creates a map plot of the given quantity. 
    figure;
    imagesc(lonT, latT, quantity)
    c = colorbar;
    ylabel(c, C_label, 'Interpreter', 'latex', 'Fontsize', Fontsize)
    set(gca, 'YDir', 'normal', 'Fontsize', Fontsize_Ticks)
    xlabel('Longitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', Fontsize)
    ylabel('Latitude ($^\circ$)', 'Interpreter', 'latex', 'Fontsize', Fontsize)
    set(gca, 'ylim', [-90 90]);
    set(gca, 'ytick', -90:30:90);
    set(gca, 'xlim', [-180 180]);
    set(gca, 'xtick', -180:30:180);
end