
addpath(genpath(pwd))

load('ReachData_3S.mat')
load('DamDatabase_3S.mat')
load('Q_3S.mat')

%% global variable initialization

%mimimum volume to be considered for mobilization
global roundpar;
roundpar = 0;
   
%% Sediment classes definition
   
sed_range = [-1.5 , 6.5]; %range of sediment sizes considered in the model
class_size = 2; %amplitude of the sediment classes

global psi
psi =  sed_range(1):class_size:sed_range(2);
clear sed_range class_size

%% Preprocessing

Network = graph_preprocessing(ReachData);

%% define timescale

timescale = 365;

%% define external sediment inputs

load('Qbi_input_3S.mat')

%% define input sediment load in the deposit layer

%initialize deposit layer
Qbi_dep_in = cell(1,length(ReachData));

for n=1:length(ReachData)
    
    Qbi_dep_in{n} = zeros(size(psi));
    
end

clear  deposit D16 D50 D84 Fi_r

%% run D-CASCADE 

[data_output,extended_output,dam_output] = DCASCADE_main_dams( ReachData , Network , Q , timescale , DamDatabase, ORparameters, Qbi_dep_in , Qbi_input , dates_Q );

%% run D-CASCADE (user defined settings)

% tr_cap_equation = 1;
% partition_formula = 3;
% velocity_formula = 2;
% 
% 
% [data_output,extended_output] = DCASCADE_main( ReachData , Network , Q , timescale , Qbi_dep_in , Qbi_input ,'tr_cap_equation',tr_cap_equation,'partition_formula',partition_formula,'velocity_formula',velocity_formula );
% 
% clear  tr_cap_equation partition_formula velocity_formula

%% plot results

dynamic_plot ( data_output, ReachData ,DamDatabase )

plot_dam_features(dam_output,Q,dates_Q)


