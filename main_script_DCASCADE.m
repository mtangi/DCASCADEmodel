% This script shows how to use the functions available in the
% MATLAB_DCASCADE repository on the case study of the Bega river system,
% NSW, Australia

%% add all folders to matlab path
addpath(genpath(pwd))

%% load input data 

% load river network as ReachData: nx1 struct reporting for each reach n of the network the attribute columun variables
load('ReachData_Bega.mat')

% load water discharge data for each network reach, for each daily timestep (m3/s)
load('Q_Bega.mat') %load discharge data for each network reach

%% Sediment classes definition
   
sed_range = [-5.5 , 3.5]; %range of sediment sizes considered in the model
class_size = 1.5; %amplitude of the sediment classes

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

%for the Bega river network example, no external sediment input is
%initialized


clear Qbi_input; Qbi_input  = cell(timescale,1); [Qbi_input{:}] = deal(zeros(length(Network.II), length(psi)));

%% define input sediment load in the deposit layer
%Qbi_dep_in is a nx1 cell structure. In each cell, it contains a 1xc matrix
%which reports the initialized sediment volume in the deposit layer for the
%first timestep of the simulation, divided into the sediment classes c 

%for the Bega river network example, the deposit in the first timestep is
%specified in the field "deposit" in reachdata [m3/km]

deposit = [ReachData.deposit].*[ReachData.Length]; %amount of material initialized as in-channel sediment deposit [m3]
% granulometry of sediment deposit for each reach
D16 = [ReachData.D16]';
D50 = [ReachData.D50]';
D84 = [ReachData.D84]';

%extract GSD from D16,D50,D84
Fi_r = GSDcurvefit( D16, D50, D84 );

%initialize deposit layer
Qbi_dep_in = cell(1,length(ReachData));

for n=1:length(ReachData)
    
    Qbi_dep_in{n} = deposit(n).*Fi_r(n,:);
    
end

clear  deposit D16 D50 D84 Fi_r

%% run D-CASCADE 

[data_output,extended_output] = DCASCADE_main( ReachData , Network , Q , timescale , Qbi_dep_in , Qbi_input  );

%% run D-CASCADE (user defined settings)
% this section shows how to run D-CASCADE with custom settings.

% tr_cap_equation = 1;
% partition_formula = 3;
% velocity_formula = 2;
% 
% [data_output,extended_output] = DCASCADE_main( ReachData , Network , Q , timescale , Qbi_dep_in , Qbi_input ,'tr_cap_equation',tr_cap_equation,'partition_formula',partition_formula,'velocity_formula',velocity_formula , 'OpRule', 3);
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
