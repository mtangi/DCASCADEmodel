function plot_time_changes(data_output , reaches , varargin)
%PLOT_TIME_CHANGES plots different subplots of the network showing the
%changes in two data in the data_plot matrix with time.
%The two dataset visualized can be chosen via the additional variable plotID 
%
% INPUT :
%
% data_output    = 10x2 Struct defining the processed output parameters of D-CASCADE
% reaches        = vector containing the ID of the n reaches which data are to be displayed (max 16 reaches)
%
%----
%
% OPTIONAL INPUT: 
% plotID       = 2x1 vector containing the ID of the two variables to be displayed, defined by their row position in data_output
% timeinterval = 2x1 vector containing the start and end time for the data displayed in the plot
% LineWidth    = value indicating the Line width, specified as a positive value
% FontSize     = value indicating the Font Size
% titles       = nx1 cell containing the title for each subplot (dafault:"Reach reachID")
% reachID      = nx1 cell containing the title for each subplot (dafault:"Reach reachID")
% showlegend   = char indicating if the legend must be show ('on' to show, 'off' to omit, default 'on')
% equal_axis   = char indicating if the two y-axis limits of each subplot must be set as equal to aid data visulization  ('on' for equal axis limits, 'off' to different axis limit, default 'off')
% colorplot    = 2x3 vector containing the RGB coloration for the two data to be displayed

%% default settings 

def_linewidth = 2;
def_cMap = 'hsv';
def_timeinterval = [];

%% read additional inputs 

p = inputParser;
addOptional(p,'plotID',[ 1 2 ]);
addOptional(p,'cMap',def_cMap);
addOptional(p,'timeinterval',def_timeinterval);
addOptional(p,'LineWidth',def_linewidth);
addOptional(p,'FontSize',12);
addOptional(p,'titles',[]);
addOptional(p,'equal_axis', 'off');
addOptional(p,'showlegend', 'on');
addOptional(p,'colorplot', [] );

parse(p,varargin{:})

plot_id = p.Results.plotID;
linewidth = p.Results.LineWidth;
TI = p.Results.timeinterval;
name_colormap = p.Results.cMap ;
Fontsize = p.Results.FontSize;
titles_subplot = p.Results.titles;
equal_axis = p.Results.equal_axis;
showlegend = p.Results.showlegend;
colorplot = p.Results.colorplot;

if isempty(TI)
    TI = [2 min(cellfun('size',data_output(plot_id,2),1))];
end
  
%% find constant colorbar for plot

cMapName = [name_colormap '(' num2str(length(data_output(:,1))) ')'];
cMap = eval(cMapName);  

if isempty(colorplot)
    colorplot = cMap(plot_id,:);
end

%% find subplot design
if length(reaches) > 12
    n_subplot = [4 4];
    reaches = reaches(1:16);
elseif length(reaches) > 9
    n_subplot = [3 4];
elseif length(reaches) > 6
    n_subplot = [3 3];
elseif length(reaches) > 4
    n_subplot = [3 2];
elseif length(reaches) > 2
    n_subplot = [2 2];
elseif length(reaches) > 1
    n_subplot = [1 2];
else
    n_subplot = [1 1];
end

%% set y axis limits, if Equal axis option is decided

if strcmp(equal_axis,'on') 
    yinterval1 = [ min( data_output{plot_id(1),2}(:,reaches ) ,[], 'all') ; max(data_output{plot_id(1),2}(:,reaches ) ,[], 'all')] ;
    yinterval2 = [ min( data_output{plot_id(2),2}(:,reaches ) ,[], 'all') ; max(data_output{plot_id(2),2}(:,reaches ) ,[], 'all')] ; 
end

%% plot data in frame
fig = figure;
set(gcf,'color','w');
set (fig, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);

for i=1:length(reaches)
    subplot(n_subplot(1),n_subplot(2),i);
   
    plot(data_output{plot_id(1),2}(:,reaches(i) ) , 'Color', colorplot(1,:) , 'LineWidth' , linewidth)

    if strcmp(equal_axis,'on') ; ylim( yinterval1 ); end
    
    yyaxis right
    hold on
    plot(data_output{plot_id(2),2}(:,reaches(i)) , 'Color', colorplot(2,:) , 'LineWidth' , linewidth)
    if strcmp(equal_axis,'on') ; ylim( yinterval2 ); end
      
    %add new title if specified
    if isempty(titles_subplot) 
        title( [ 'Reach ' num2str(reaches(i))] );
    else      
        title( [titles_subplot{i} ] );
    end

    hold off 
    xlim(TI)

    grid minor

    set(gca,'XGrid','on','FontSize',Fontsize)
end

%% set legend

if strcmp(showlegend,'on')
    legend( data_output(plot_id,1),'FontSize',15,'location', 'southeast')
end

end
