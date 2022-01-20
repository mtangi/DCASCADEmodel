function [] = plot_dam_features(dam_output,Q,dates_Q)


%% load dam features

figure
set(gcf,'color','w');
DamDatabase = dam_output{1,2};
ResVolume = dam_output{3,2};
release = dam_output{4,2};


n_dams = length(DamDatabase);

ts = size(ResVolume,1)-2;

%% plot dam features

dam_heigth = zeros(ts,size(release,2));

for d=1:n_dams
    sb = subplot(n_dams,1,d);
    pos = get(sb, 'position');

    plot( datetime(dates_Q(:,1:ts)'),[Q(1:ts, DamDatabase(d).Node_ID)])
    hold on

    %dam name annotation
    str = ([ [DamDatabase(d).Code]]); 
    
    annotation('textbox',pos,'String',str,'FitBoxToText','on','EdgeColor', 'none','FontWeight','bold','FontSize',14)

    %plot release (outflow)
    plot( datetime(dates_Q(:,1:ts)'),[release(1:ts,d)])
    hold on

    spillway = max(0,release(1:ts,d) - DamDatabase(d).design_discharge);
    plot( datetime(dates_Q(:,1:ts)'),spillway)

    ylim([0, max([release;Q(1:ts, [DamDatabase.Node_ID])] ,[],'all')]);
    ylabel('Discharge [m^3/s]')

    yyaxis right

    dam_heigth(:,d) = Reservoir_LSWConversion( [ResVolume(1:ts,d)] ,DamDatabase(d).WL_table, 3, 1);
 
    ylim([0, max(dam_heigth,[],'all')+0.5]); % for the 1 dam simulation


    plot( datetime(dates_Q(:,1:ts)'),[dam_heigth(:,d)],'LineWidth',1.5)
    %plot(datetime(dates_Q(:,1:ts)'),[ResVolume(1:ts,d)])

    %title([DamDatabase_active(d).Name ' - reach ' num2str(DamDatabase_active(d).Node_ID)])

    ylabel('Water level [m]')
    set(gca,'FontSize',14)

    grid on

    set(gca, 'YGrid', 'off', 'XGrid', 'on')
end

legend({'Inflow [m^3/s]','Outflow [m^3/s]','Spillway [m^3/s]','Reservoir heigth [m]'},'Location','southoutside','Orientation','horizontal')
set(gcf,'color','w');

end

