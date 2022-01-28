function [] = dynamic_plot_Q ( data_output, ReachData , varargin )
%plot output data of the D-CASCADE model 

%% read additional inputs 

Q = data_output{cell2mat(cellfun(@(x)strcmp(x,'Discharge'),data_output(:,1), 'UniformOutput', false)) == 1 ,2} ;

def_ShowID = 'off';

p = inputParser;
addOptional(p,'dates_Q',size(Q,1));

parse(p,varargin{:})

dates_Q = p.Results.dates_Q ;

%% input data

start_time = 2;
plot_class = 3; %default plot variable
cMap = 'jet';

sim_length = min(cell2mat(cellfun(@(x)size(x,1),data_output(1:length(data_output) - 1,2), 'UniformOutput', false)));

%% define plot variables

%define sediment classes
n_class = 30;
i_class = 100/n_class - 0.00001; %interval between classestext

for c=1:length(data_output) - 1

    cClass{c} = unique(prctile(data_output{c,2}(data_output{c,2}~=0),0:i_class:100)); 
%     if c == length(plot_var)
%          cClass{c} = 0:5:100;
%     end

end

if any(cell2mat(cellfun(@(x)strcmp(x,'Total sed in the reach - per class'),data_output(:,1), 'UniformOutput', false)))
   % cClass{ cell2mat(cellfun(@(x)strcmp(x,'SedDelRt'),data_output(:,1), 'UniformOutput', false)) == 1 } = 0:10:100;
end

%% plot starting river network

figure
subplot(4,1,[1,2,3])

h =  findobj('type','figure'); %find open figures
n = length(h);
fig = figure(n); %enlarge figure (the last opened)
set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);

plotvariable = data_output{plot_class,2}(start_time,:);

plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});

%plot discharge
outlet = find([ReachData.FromN] == [ReachData.ToN]); %I show the discharge referred to the outlet
qwindow = 120; %range of discharge values to be visualized

subplot(4,1,[4])



qp = plot( datetime(dates_Q(:,1:size(Q,1))'),Q(: , outlet),'b') ;
hold on 
dl = plot([ datetime(dates_Q(:,1)') datetime(dates_Q(:,1)') ], [0 100000],'r','LineWidth',3);
ylim([0, max(Q(:, outlet))+100]);
xlim([datetime(dates_Q(:,1)') - days(qwindow), datetime(dates_Q(:,1)') + days(qwindow)]);
ylabel('Discharge [m^3]')

%% display reach ID, nodes , dams and additional sed fluxes

%display reach nodes
outlet_ID = find([ReachData.FromN] == [ReachData.ToN]);

hold on
hn = scatter([ReachData.x_FN ReachData(outlet_ID).x_TN],  [ReachData.y_FN ReachData(outlet_ID).y_TN] ,30,'o','filled','k');
set(hn,'DisplayName','Nodes','Visible','off','HandleVisibility','off');

%display reach ID
for i=1:size(ReachData,1)
    str{i} = num2str([ReachData(i).reach_id]);
end

xt = ([ReachData.x_FN]+[ReachData.x_TN])/2;
yt = ([ReachData.y_FN]+[ReachData.y_TN])/2;
hold on
ts = textscatter(xt,yt,str);
set(ts,'Visible','off','MarkerColor','none', 'TextDensityPercentage' ,80 ,'HandleVisibility','off');

%% button code initialitation

default_button = double('0');
node_button = double('n');
ID_button = double('i');
int_button = double('1');
sel_button = double('2');
list_button = double('3');
for_button = 29;
back_button = 28;

def_plot_class = 1;

timestep = start_time;
timelapse = 10;

exit_button = [27 8]; % exitbutton = ESC or Backsapce

option_button = default_button;

% str = {[ char(int_button) ' : reach selection via network'],...
%     [char(man_button) ' : reach selection via reach ID'],[char(ID_button) ' : show reach ID'], ...
%     [char(node_button) ' : show network nodes']};
% 
% str{length(str)+1} = 'ESC or BS : close figure';
% 
% annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
%     'BackgroundColor', 'w','FontSize',12);
% 

str = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
 'BackgroundColor', 'w','FontSize',12);
 
%% plot detail on the network reach
while ~ any (option_button == exit_button)  
    
    % load use command
    w = waitforbuttonpress;

    if w

        p = get(gcf, 'CurrentCharacter');
        option_button = double(p);
        %displays the ascii value of the character that was pressed    
    end

    %go forward in time
    if option_button == for_button

      subplot(4,1,[1,2,3])

      timestep = min(timestep+timelapse, sim_length );

      old = findall(gcf,'Type','line');
      delete(old);

      plotvariable = data_output{plot_class,2}(timestep,:);
      plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});

      delete(a)
      str = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
      a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
     'BackgroundColor', 'w','FontSize',12);

      uistack(ts,'top'); uistack(hn,'top');      %put nodes and reach ID on top of the new plotted network

      subplot(4,1,4)

      delete(dl)
      delete(qp)
      qp = plot( datetime(dates_Q'),Q(:, outlet),'b') ;
      hold on 
      dl = plot([ datetime(dates_Q(:,timestep)') datetime(dates_Q(:,timestep)')], [0 100000],'r','LineWidth',3);

      xlim([datetime(dates_Q(:,timestep)') - days(qwindow), datetime(dates_Q(:,timestep)') + days(qwindow)]);

    %go backward in time
    elseif option_button == back_button

      subplot(4,1,[1,2,3])

      timestep = max(timestep-timelapse, 1);

      old = findall(gcf,'Type','line');
      delete(old);

      plotvariable = data_output{plot_class,2}(timestep,:);
      plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});

      delete(a)
      str = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
      a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
     'BackgroundColor', 'w','FontSize',12);

      uistack(ts,'top'); uistack(hn,'top');      %put nodes and reach ID on top of the new plotted network

      subplot(4,1,4)

      delete(dl)
      delete(qp)
      qp = plot( datetime(dates_Q'),Q(:, outlet),'b') ;
      hold on 
      dl = plot([ datetime(dates_Q(:,timestep)') datetime(dates_Q(:,timestep)')], [0 100000],'r','LineWidth',3);

      xlim([datetime(dates_Q(:,timestep)') - days(qwindow), datetime(dates_Q(:,timestep)') + days(qwindow)]);

    %change timelapse
    elseif option_button == int_button

      answer = inputdlg('Enter timelapse','Input timelapse',[1 35],{num2str(timelapse)});
      if ~isempty(answer)
         timelapse = str2double(answer{1});
      end

      delete(a)
      str = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
      a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
     'BackgroundColor', 'w','FontSize',12);

    %change timestep
    elseif option_button == sel_button

        answer = inputdlg('Enter timestep','Input timestep',[1 35],{num2str(timestep)});

        if ~isempty(answer)

         timestep = min(str2double(answer{1}), sim_length);
        end

        subplot(4,1,[1,2,3])

        old = findall(gcf,'Type','line');
        delete(old);

        plotvariable = data_output{plot_class,2}(timestep,:);
        plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});

        delete(a)
        str = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' 1] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
        a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
        'BackgroundColor', 'w','FontSize',12);

        subplot(4,1,4)

        delete(dl)
        delete(qp)
        qp = plot( datetime(dates_Q'),Q(:, outlet),'b') ;
        hold on 
        dl = plot([ datetime(dates_Q(:,timestep)') datetime(dates_Q(:,timestep)')], [0 100000],'r','LineWidth',3);

        xlim([datetime(dates_Q(:,timestep)') - days(qwindow), datetime(dates_Q(:,timestep)') + days(qwindow)]);

    %change plot_variable
    elseif option_button == list_button

        %old_plot_class = plot_class;
        %ask user to select plot_variable

        [plot_class,tf] = listdlg('ListString',data_output(:,1),'PromptString','Select plot_variable:',...
                            'SelectionMode','single','InitialValue',1,'ListSize',[250,250],'CancelString','Default');
        if tf == 0
          plot_class = def_plot_class;
        end

        plotvariable = data_output{plot_class,2}(timestep,:);

        % plot new plot variable
        old = findall(gcf,'Type','line');
        delete(old);

        plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});


        %change annotation
        delete(a)
        str = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' 1] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
        a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str,'FitBoxToText','on',...
        'BackgroundColor', 'w','FontSize',12);

        uistack(ts,'top'); uistack(hn,'top');      %put nodes and reach ID on top of the new plotted network

    %toggle nodes visibility
    elseif option_button == node_button 

       if strcmp(get(hn,'visible'), 'on')
             set(hn,'Visible','off','HandleVisibility','off')
       else
           set(hn,'Visible','on','HandleVisibility','on')
       end

    %toggle additional sed flows visibility
    elseif nargin > 10 &&  ~isempty(additional_sed_flow) && option_button == add_button 

       if strcmp(get(ha,'visible'), 'on')
             set(ha,'Visible','off')
       else
           set(ha,'Visible','on')
       end

    %toggle node id
    elseif option_button == ID_button 

       if strcmp(get(ts,'visible'), 'on')
             set(ts,'Visible','off')
       else
           set(ts,'Visible','on')
           uistack(ts,'top');

       end

    end 
   
   
end

close

end