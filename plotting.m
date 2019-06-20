function plotting(path, figpath, file, filter, dlist, hlist, trans)
    return % to plot nothing
    if filter == 'H'
        M = dlmread(strcat(figpath,sprintf('%s-spectrum.txt', filter)));
        lambda = M(:,1);
        spec = M(:,2);
        sigma = M(:,3);
        sspec = M(:,4);
        sflat = M(:,5);
        wavenumber = M(:,6);
        nc_spec31 = M(:,7);
        nc_spec42 = M(:,8);
        nc_spec53 = M(:,9);
        %obsTime = M(:,6); % not ready yet
        
        Q_31_range = 105:161;
        Q_42_range = 350:426; % Q branch with linear noise
        Q_53_range = 618:698;
        
        figure;
        plot(lambda(Q_31_range), nc_spec31(Q_31_range), 'b-', lambda(Q_31_range), lambda(Q_31_range)*0, 'r:')
        title(sprintf('%s-%s Noise Cancelled at 31 Q-branch', filter, file));
        xlabel('\lambda [{\mu}m]');
        ylabel('Integrated intensity [ADU]');
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambda(114) lambda(114)],y1, '--m', 'linewidth', 1)
        plot([lambda(133) lambda(133)],y1, '--m', 'linewidth', 1)
        hold off
          
        figure;
        plot(lambda(Q_42_range), nc_spec42(Q_42_range), 'b-', lambda(Q_42_range), lambda(Q_42_range)*0, 'r:')
        title(sprintf('%s-%s Noise Cancelled at 42 Q-branch', filter, file));
        xlabel('\lambda [{\mu}m]');
        ylabel('Integrated intensity [ADU]');
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambda(381) lambda(381)],y1, '--m', 'linewidth', 1)
        plot([lambda(400) lambda(400)],y1, '--m', 'linewidth', 1)
        hold off
        
        figure;
        plot(lambda(Q_53_range), nc_spec53(Q_53_range), 'b-', lambda(Q_53_range), lambda(Q_53_range)*0, 'r:')
        title(sprintf('%s-%s Noise Cancelled at 53 Q-branch', filter, file));
        xlabel('\lambda [{\mu}m]');
        ylabel('Integrated intensity [ADU]');
        y1=get(gca,'ylim');
        hold on % Vertical lines to visualize integration area
        plot([lambda(652) lambda(652)],y1, '--m', 'linewidth', 1)
        plot([lambda(672) lambda(672)],y1, '--m', 'linewidth', 1)
        hold off
    else
        M = dlmread(strcat(figpath,sprintf('%s-spectrum.txt', filter)));
        lambda = M(:,1);
        spec = M(:,2);
        sigma = M(:,3);
        sspec = M(:,4);
        sflat = M(:,5);
        wavenumber = M(:,6);
        %obsTime = M(:,6); % not ready yet 
    end % if filter == 'H'
    
    figure;
    plot(lambda,spec);
    title(sprintf('%s-%s spectrum', filter, file));
    xlabel('\lambda [{\mu}m]');
    ylabel('Integrated intensity [ADU]');
%     
    figure;
    plot(wavenumber, spec)
    title(sprintf('%s-%s spectrum', filter, file));
    xlabel('Wavenumber [cm^{-1}]');
    ylabel('Integrated intensity [ADU]');
   

    % Plotting for debugging/learning

%     M = dlmread(strcat(figpath,sprintf('%s-spectrum.txt', filter)));
%     lambda=M(:,1);
%     spec=M(:,2);
%     spectrum = figure;
%     plot(lambda, spec);
%     title("Spectrum")
%     xlabel("Wavelength [\mu m]")
%     ylabel("Intensity")
%     %print(strcat(figpath, '\\spectrum.pdf'), '-dpdf', '-bestfit')
% 
%     data_unproc_sat1000 = dlmread(strcat(figpath, 'data_unproc_sat1000.txt'));
%     data_unproc = figure;
%     contourf(data_unproc_sat1000,'LineColor','none');
%     pcolor(data_unproc_sat1000); shading interp;
%     colormap('hot(600)')
%     colorbar;
%     title("data, unprocessed")
%     xlabel('Point on the sky')
%     ylabel('Wavelength')
%     %print(strcat(figpath, '\\unproc_sat150.pdf'), '-dpdf', '-bestfit')
% 
%     data_plot = dlmread(strcat(figpath, 'data_preflat_sat150.txt'));
%     data_preflat = figure;
%     contourf(data_plot,'LineColor','none');
%     pcolor(data_plot); shading interp;
%     colormap('hot(600)')
%     colorbar;
%     title("data, filtered")
%     xlabel('Point on the sky')
%     ylabel('Wavelength')
% 
%     %print(strcat(figpath, '\\preflat_sat150.pdf'), '-dpdf', '-bestfit')
% 
%     data_plot = dlmread(strcat(figpath, 'data_postflat_sat150.txt'));
%     data_postflat = figure;
%     contourf(data_plot,'LineColor','none');
%     pcolor(data_plot); shading interp;
%     colormap('hot(600)')
%     colorbar;
%     title("data, filtered")
%     xlabel('Point on the sky')
%     ylabel('Wavelength')
% 
%     %print(strcat(figpath, '\\postflat_sat150.pdf'), '-dpdf', '-bestfit')
% 
%     data_plot = dlmread(strcat(figpath, 'data_lin_sat150.txt'));
%     data_linearized = figure;
%     contourf(data_plot,'LineColor','none');
%     pcolor(data_plot); shading interp;
%     colormap('hot(600)')
%     colorbar;
%     title("data, linearized")
%     xlabel('Point on the sky')
%     ylabel('Wavelength')
%     
    %print(strcat(figpath, '\\linearized_sat150.pdf'), '-dpdf', '-bestfit')
% 
% %   Save figures
% 
% %     savefig(data_unproc, sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s_figs\\dataUnproc.fig',file))
% %     savefig(data_preflat, sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s_figs\\dataPreflat.fig',file))
% %     savefig(data_postflat, sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s_figs\\dataPostflat.fig',file))
% %     savefig(data_linearized, sprintf('C:\\Users\\birke\\Documents\\prosjekt\\ReduceSingleNotfileBirk\\%s_figs\\dataLinearized.fig',file))

end

    