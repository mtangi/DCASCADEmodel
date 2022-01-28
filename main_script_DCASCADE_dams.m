
% This script shows how to use the functions available in the
% MATLAB_DCASCADE repository on the case study of the 3S river system,
% a tributary of the Mekong.
%
% This case study include multiple reservoirs in the system, whose features
% (stored water and sediment volumes, release, flooded area and others) are
% dynamically modelle alongsize sediment delivery and transport 

%% add all folders to matlab path
addpath(genpath(pwd))

%% load input data 

% load river network as ReachData: nx1 struct reporting for each reach n of the network the attribute column variables
load('ReachData_3S.mat')

% load water discharge data for each network reach, for each daily timestep (m3/s)
load('Q_3S.mat') %load discharge data for each network reach

% load dam information as DamDatabase: nx1 struct reporting for each reservoir d its location and features 
load('DamDatabase_3S.mat')

%% Sediment classes definition
   
sed_range = [-1.5 , 6.5]; %range of sediment sizes considered in the model
class_size = 2; %amplitude of the sediment classes

global psi
psi =  sed_range(1):class_size:sed_range(2);
clear sed_range class_size

%% Preprocessing
% graph_preprocessing extract network connectivity informations from the network structure 
Network = graph_preprocessing(ReachData);

%% define timescale 
% timescale represent the length of simulation

timescale = 365;

%% define external sediment inputs
%Qbi_input is a tx1 cell structure. In each cell, it contains a nxc matrix
%which reports the input sediment volume in timestep t from external
%sources (e.g.  hillslope erosion or debris flows) for each reach n,
%divided into the sediment classes c 

%for the 3S river network example, the external sediment inputs are already
%defined based on sediment yield informations from Kondolf et al., 2014

load('Qbi_input_3S.mat')

%% define input sediment load in the deposit layer
%Qbi_dep_in is a nx1 cell structure. In each cell, it contains a 1xc matrix
%which reports the initialized sediment volume in the deposit layer for the
%first timestep of the simulation, divided into the sediment classes c 

%for the 3S river network example, no deposit is initialized.

%initialize deposit layer
Qbi_dep_in = cell(1,length(ReachData));

for n=1:length(ReachData)
    
    Qbi_dep_in{n} = zeros(size(psi));
    
end

clear  deposit D16 D50 D84 Fi_r

%% run D-CASCADE 

[data_output,extended_output,dam_output] = DCASCADE_main_dams( ReachData , Network , Q , timescale , DamDatabase, ORparameters, Qbi_dep_in , Qbi_input ,'dates_Q', dates_Q );

%% run D-CASCADE (user defined settings)
% this section shows how to run D-CASCADE with custom settings.
% 
% tr_cap_equation = 1;
% partition_formula = 3;
% velocity_formula = 2;
% operating_rule = 3; %specify operating rule
% 
% [data_output,extended_output, dam_output] = DCASCADE_main( ReachData , Network , Q , timescale , Qbi_dep_in , Qbi_input , ...
%     'OpRule', operating_rule,'tr_cap_equation',tr_cap_equation,'partition_formula',partition_formula,'velocity_formula',velocity_formula);
% 
% clear  tr_cap_equation partition_formula velocity_formula

%% plot results

%dynamic plot displays the results collected in data output.
%this interactive function allows the user to move forward or backward in
%time, plot different outputs and explore how these output changes over
%time
dynamic_plot ( data_output, ReachData  )

% plot_time_changes shows the evolution through time of two parameters in
% data_output (specified as input plotID, where the ID is the position of
% the parameter in data_output ), for a number of user-specified reaches,
% defined as second input

plot_time_changes(data_output, 94:99 , 'plotID', [3 9] )

% plot_dam_features shows the evolution through time of the reservoir
% features, including water level, water input, standard release and
% release through spillways

plot_dam_features(dam_output,Q,dates_Q)
