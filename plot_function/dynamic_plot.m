function [] = dynamic_plot ( data_output, ReachData , DamDatabase )
%DYNAMIC_PLOT displays the results collected in data output.
%This interactive function allows the user to move forward or backward in
%time, plot different outputs and explore how these output changes over
%time.

%% input data

global psi
dmi = 2.^(-psi);

start_time = 2;
plot_class = 3; %default plot variable
cMap = 'jet';

sim_length = min(cell2mat(cellfun(@(x)size(x,1),data_output(1:length(data_output) - 1,2), 'UniformOutput', false)));

if nargin <3  
    DamDatabase_active = [];
else
    DamDatabase_active = DamDatabase([DamDatabase.portfolio] == 1) ;
end

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

%extract tot_sed_class values
if any(cell2mat(cellfun(@(x)strcmp(x,'Total sed in the reach - per class [m^3]'),data_output(:,1), 'UniformOutput', false)))
    
    def_sed_class = 1;
    
    % find position of tot_sed_class in data_output
    pos_tot_sed_class = find(cell2mat(cellfun(@(x)strcmp(x,'Total sed in the reach - per class [m^3]'),data_output(:,1), 'UniformOutput', false)) ==1);
    %load tot_sed_class
    tot_sed_class = data_output{pos_tot_sed_class,2};

    %ubstitute tot_sed_class with tot_sed_class for the def_class in data_output
    data_output{pos_tot_sed_class,2} = tot_sed_class{def_sed_class};
    data_output{pos_tot_sed_class,1} = ['Total sed - per class '];
    
    cClass{pos_tot_sed_class} = unique(prctile(data_output{pos_tot_sed_class,2}(data_output{pos_tot_sed_class,2}~=0),0:i_class:100)); 

else
    pos_tot_sed_class = [];
end

Q = data_output{cell2mat(cellfun(@(x)strcmp(x,'Discharge [m^3/s]'),data_output(:,1), 'UniformOutput', false)) == 1 ,2} ;

%% plot starting river network

figure
h =  findobj('type','figure'); %find open figures
n = length(h);
fig = figure(n); %enlarge figure (the last opened)
set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);

plotvariable = data_output{plot_class,2}(start_time,:);

plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});

%% display reach ID, nodes , dams and additional sed fluxes

%display reach nodes
outlet_ID = find([ReachData.FromN] == [ReachData.ToN]);

hold on
hn = scatter([ReachData.x_FN ReachData(outlet_ID).x_TN],  [ReachData.y_FN ReachData(outlet_ID).y_TN] ,30,'o','filled','k');
set(hn,'DisplayName','Nodes','Visible','off','HandleVisibility','off');

%display reach ID
for i=1:size(ReachData,1)
    str_a{i} = num2str([ReachData(i).reach_id]);
end

xt = ([ReachData.x_FN]+[ReachData.x_TN])/2;
yt = ([ReachData.y_FN]+[ReachData.y_TN])/2;
hold on
ts = textscatter(xt,yt,str_a);
set(ts,'Visible','off','MarkerColor','none', 'TextDensityPercentage' ,80 ,'HandleVisibility','off');

%display Dams names
if ~isempty(DamDatabase_active)

    xt = ([ReachData([DamDatabase_active.Node_ID]).x_FN]);
    yt = ([ReachData([DamDatabase_active.Node_ID]).y_FN]);
    hold on
    dm = scatter(xt,yt,300,'v','filled','r');
    set(dm,'DisplayName','Dams','Visible','on','HandleVisibility','off');
    
else
    
    dm = scatter([ReachData(1).x_FN],[ReachData(1).y_FN]);
    set(dm,'Visible','off','MarkerFaceColor','none','HandleVisibility','off');
    
end
%% button code initialitation

default_button = double('0');
node_button = double('n');
dam_button = double('d');
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

str_a = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['variable = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str_a,'FitBoxToText','on',...
 'BackgroundColor', 'w','FontSize',12);
 

str_b = {['right arrow - go forward in time' ] ['left arrow - go backward in time' ] ['1 - input timelapse' ] ['2 - input timestep'] ['3 - choose plot variable'] ['i -show reach ID'] ['n -show network nodes']};
b = annotation('textbox',[0.12 0.62 0.3 0.3],'String',str_b,'FitBoxToText','on','BackgroundColor', 'w','FontSize',12);

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

      timestep = min(timestep+timelapse, sim_length );
      
      old = findall(gcf,'Type','line');
      delete(old);
      
      plotvariable = data_output{plot_class,2}(timestep,:);
      plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});
       
      delete(a)
      str_a = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
      a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str_a,'FitBoxToText','on',...
     'BackgroundColor', 'w','FontSize',12);
      b = annotation('textbox',[0.12 0.62 0.3 0.3],'String',str_b,'FitBoxToText','on','BackgroundColor', 'w','FontSize',12);

      uistack(ts,'top'); uistack(hn,'top'); uistack(dm,'top');     %put nodes and reach ID on top of the new plotted network
      
   %go backward in time
   elseif option_button == back_button
       
      timestep = max(timestep-timelapse, 1);
      
      old = findall(gcf,'Type','line');
      delete(old);
      
      plotvariable = data_output{plot_class,2}(timestep,:);
      plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});
        
      delete(a)
      str_a = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
      a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str_a,'FitBoxToText','on',...
     'BackgroundColor', 'w','FontSize',12);
      b = annotation('textbox',[0.12 0.62 0.3 0.3],'String',str_b,'FitBoxToText','on','BackgroundColor', 'w','FontSize',12);

      uistack(ts,'top'); uistack(hn,'top'); uistack(dm,'top');     %put nodes and reach ID on top of the new plotted network
      
   %change timelapse
   elseif option_button == int_button
       
      answer = inputdlg('Enter timelapse','Input timelapse',[1 35],{num2str(timelapse)});
      if ~isempty(answer)
         timelapse = str2double(answer{1});
      end
     
      delete(a)
      str_a = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
      a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str_a,'FitBoxToText','on',...
     'BackgroundColor', 'w','FontSize',12);
      b = annotation('textbox',[0.12 0.62 0.3 0.3],'String',str_b,'FitBoxToText','on','BackgroundColor', 'w','FontSize',12);
   
    %change timestep
    elseif option_button == sel_button

        answer = inputdlg('Enter timestep','Input timestep',[1 35],{num2str(timestep)});

        if ~isempty(answer)

         timestep = min(str2double(answer{1}), sim_length);
        end

        old = findall(gcf,'Type','line');
        delete(old);

        plotvariable = data_output{plot_class,2}(timestep,:);
        plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});

        delete(a)
        str_a = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' 1] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
        a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str_a,'FitBoxToText','on',...
        'BackgroundColor', 'w','FontSize',12);
        b = annotation('textbox',[0.12 0.62 0.3 0.3],'String',str_b,'FitBoxToText','on','BackgroundColor', 'w','FontSize',12);
    
        uistack(ts,'top'); uistack(hn,'top'); uistack(dm,'top');     %put nodes and reach ID on top of the new plotted network

    %change plot_variable
    elseif option_button == list_button
       
        %old_plot_class = plot_class;
        %ask user to select plot_variable
        data_output{pos_tot_sed_class,1} = ['Total sed - per class '];

        [plot_class,tf] = listdlg('ListString',data_output(:,1),'PromptString','Select plot_variable:',...
                            'SelectionMode','single','InitialValue',1,'ListSize',[250,250],'CancelString','Default');
        if tf == 0
          plot_class = def_plot_class;
        end
        
        % if I select tot_sed_class, I choose wich class to display
        if plot_class == pos_tot_sed_class
            [sed_class,tf] = listdlg('ListString',string(dmi),'PromptString',' Select grain size class [mm]:',...
                        'SelectionMode','single','InitialValue',1,'ListSize',[250,250],'CancelString','Default');
            if tf == 0
                sed_class = def_sed_class;
            end          
            
            data_output{pos_tot_sed_class,2} = tot_sed_class{sed_class};
            data_output{pos_tot_sed_class,1} = ['Total sed in the reach - class ' , num2str(dmi(sed_class))];
            
            cClass{pos_tot_sed_class} = unique(prctile(data_output{pos_tot_sed_class,2}(data_output{pos_tot_sed_class,2}~=0),0:i_class:100)); 

        end

        plotvariable = data_output{plot_class,2}(timestep,:);

        % plot new plot variable
        old = findall(gcf,'Type','line');
        delete(old);

        plot_network_dyn ( ReachData, plotvariable, 'CMap', cMap, 'cClass',cClass{plot_class});


        %change annotation
        delete(a)
        str_a = {['time = ' num2str(timestep)] [] ['timelapse = ' num2str(timelapse)] ['class = ' data_output{plot_class,1}] [] [' Q = ' num2str(max(Q(timestep,:)))]   };
        a = annotation('textbox',[0.841 0.62 0.3 0.3],'String',str_a,'FitBoxToText','on',...
        'BackgroundColor', 'w','FontSize',12);
        b = annotation('textbox',[0.12 0.62 0.3 0.3],'String',str_b,'FitBoxToText','on','BackgroundColor', 'w','FontSize',12);
    
        uistack(ts,'top'); uistack(hn,'top'); uistack(dm,'top');     %put nodes and reach ID on top of the new plotted network
        
    %toggle nodes visibility
    elseif option_button == node_button 

       if strcmp(get(hn,'visible'), 'on')
             set(hn,'Visible','off','HandleVisibility','off')
       else
           set(hn,'Visible','on','HandleVisibility','on')
       end
       
   %toggle reach id
   elseif option_button == ID_button 
      
       if strcmp(get(ts,'visible'), 'on')
             set(ts,'Visible','off')
       else
           set(ts,'Visible','on')
           uistack(ts,'top');

       end
       
   %toggle reach id
   elseif and(option_button == dam_button , ~isempty(DamDatabase_active))
      
       if strcmp(get(dm,'visible'), 'on')
             set(dm,'Visible','off')
       else
           set(dm,'Visible','on')
           uistack(dm,'top');

       end
       
    
   end 
   
   

end

close

end

