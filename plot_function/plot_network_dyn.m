function plot_network_dyn(ReachData,  plotvariable, varargin)
% PLOT_NETWORK_DYN plots the river network and visualizes CONTINOUS DATA.  
%
% Can use color to represent different attributes. 
% Colorcode can be clarified using a legend. 
% Line width can represent different attributes.

% INPUT:
%
% ReachData: Nx1 struct containing from- and to-node infomation and attribute values. 
% plotvariable : Nx1 vector containing for N reaches the values to be displayed.
%
% Optional input (name-value pair argument):
%
% ShowID.: if 'on', displays the reach ID on the midpoint of the reach
% cMap: name of color map used in the figure (default : parula)
% legend_type : if 'colorbar', displays data with colorbar, else
%                   the data will be displayed with percentiles color classes
% title : Name of the output figure;
% Linewidth : width of the reach lines, either a single value of a Nx1
%                   vector of reach attributes;
% ClassNumber: number of percentiles color classes 
% cClass : vector of values of classes used for the legend

%----
%
% OUTPUT
% 
% fig = Figure object. Use f to query or modify properties of the figure after it is created.
%

%% default settings 

def_ClassNumber = 10;
def_linewidth = 4;
def_cMap = 'parula';
def_ShowID = 'off';
def_legendtype = 'percentile';

fontsize = 12;
title_fontsize = 12;

%function needs 1xN vector, switch to 1xN if vector is Nx1
if nargin ~= 1 && isequal(size(plotvariable), [length(ReachData), 1])
    plotvariable = plotvariable';
end

if nargin == 1
    plotvariable = [];
end

%% read additional inputs 

p = inputParser;
addOptional(p,'ShowID',def_ShowID);
addOptional(p,'cMap',def_cMap);
addOptional(p,'title',[]);
addOptional(p,'legendon',1);
addOptional(p,'legendtype',def_legendtype);
addOptional(p,'LineWidth',def_linewidth);
addOptional(p,'ClassNumber',def_ClassNumber);
addOptional(p,'cClass',[]);
addOptional(p,'plotnodes',[]);

parse(p,varargin{:})

indx_showID = p.Results.ShowID ;
name_colormap = p.Results.cMap ;
figname = p.Results.title ;
legendtype = p.Results.legendtype ;
n_class = p.Results.ClassNumber;
legendon = p.Results.legendon;

%% set width, color Class and plot nodes

%set width value for all reaches
if length(p.Results.LineWidth) ~= length(ReachData)
   line_width = repmat(p.Results.LineWidth(1),length(ReachData),1);
else
   line_width = p.Results.LineWidth;
end


if isempty(p.Results.cClass)
    i_class = 100/n_class + 0.001; %interval between classestext
    cClass = unique(prctile(plotvariable(plotvariable~=0),0:i_class:100)); 
    if ~isempty(plotvariable); cClass(end) = max(plotvariable); end
    cClass = [0, cClass(cClass~=0)]; %leave just 1 class equal to 0
  
else
   cClass = p.Results.cClass;
end

% find network to be plotted give the node ID as input
if ~isempty(p.Results.plotnodes)
    plotvariable = plotvariable(p.Results.plotnodes);
    ReachData = ReachData(p.Results.plotnodes);
end
  

%% plot "empty "network
          
if nargin == 1 || isempty(plotvariable)

    for i=1:length(ReachData)       % plot entire river network  
        plot([ReachData(i).X],[ReachData(i).Y],'Color',[0.7 0.7 0.7],'LineWidth',line_width(i))
        hold on
    end
     
%% plot network with data and colorbar

elseif strcmp(legendtype,'colorbar')
    
    %read colormap
    cMapName=[name_colormap '( 100 )'];
    cMap=eval(cMapName);  
    
    %find color in colormap
    pos = ceil((plotvariable-min(plotvariable))./(max(plotvariable)-min(plotvariable)).*99+0.001);

    tick = 0:0.1:1;
    tick_val = tick.*(max(plotvariable)-min(plotvariable))+min(plotvariable);
    
    figure;
    
    for i=1:length(ReachData)
        
        plot([ReachData(i).X],[ReachData(i).Y],'Color',cMap(pos(i),:),'LineWidth',line_width(i));
        hold on
        
    end
    
    % define the visualization of plot values for the legend
    if max(plotvariable)>10^4; legEnt=num2str(tick_val','%10.2E'); else legEnt=num2str(tick_val','%.4g'); end; 

    
    colormap(name_colormap)
    colorbar('Ticks',tick,...
             'TickLabels',legEnt);
            
%% plot network with data and color classes

else
          
    %define color map
    cMapLength = length(unique(cClass))+1;
    
    %read colormap
    cMapName = [name_colormap '(' num2str(cMapLength) ')'];
    cMap = eval(cMapName);  
    
  % loop through all classes 
    for c_cl=1:cMapLength
     
      % find all observations that have an attribute value that falls
      % within the current c class
      if c_cl==1
          cClassMem = find(plotvariable<=cClass(c_cl)); 
      elseif c_cl<=length(cClass) && c_cl>1
           cClassMem = find(plotvariable>cClass(c_cl-1) & plotvariable<=cClass(c_cl));
      else
          cClassMem = find(plotvariable>cClass(c_cl-1));
      end 
%       catt(catt==0)=nan; 
          
        for ll=(cClassMem) %=find(raw_data(:,ID_FromN)>0)'; % plot entire river network
        
        plot([ReachData(ll).X],[ReachData(ll).Y],'Color',cMap(c_cl,:),'LineWidth',line_width(ll));
        hold on
        
        end
    
    end
    
   %create fake lines for the legend. 
   for leggg=1:length(unique(cClass))
        hh_leg(leggg) = line([ReachData(1).x_FN ReachData(1).x_FN],[ReachData(1).y_FN ReachData(1).y_FN],'color',cMap(leggg,:),'linewidth',3);  % for each color attribute value make a fake line with length 0. 
   end
   
   if max(cClass)>10^4; legEnt=num2str(cClass','%10.2E'); else legEnt=num2str(cClass','%.4g'); end; % prepare the text for the legend
    
   % create legend 
   if legendon == 1
       leg = legend(hh_leg,legEnt);
       set(leg,'Location', 'southeastoutside');
       set(leg,'FontSize',12);
       leg_pos=get(leg,'Position');
       leg.FontSize = fontsize;
   end
    
  
end

%% enlarge figure


set(gca,'xtick',[])
set(gca,'ytick',[])
% 
% h =  findobj('type','figure'); %find open figures
% n = length(h);
% fig = figure(n); %enlarge figure (the last opened)
% set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);

%set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.5]);

%% add axes name and change font size

%xlabel('X coordinates');
%ylabel('Y coordinates');
pax = gca;
pax.FontSize = fontsize;
pax.FontWeight = 'bold';

t = title(figname);
t.FontSize = title_fontsize;

   %% plot node ID

if strcmp(indx_showID,'on')
    
    for i=1:size(ReachData,1)
        str{i} = num2str([ReachData(i).reach_id]);
    end
    xt = ([ReachData.x_FN]+[ReachData.x_TN])/2;
    yt = ([ReachData.y_FN]+[ReachData.y_TN])/2;
    
    t = textscatter(xt,yt,str, 'MarkerColor','none');
    set(t,'Visible','on','MarkerColor','none', 'TextDensityPercentage' ,60 ,'HandleVisibility','off');

end

