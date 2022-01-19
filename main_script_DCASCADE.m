
addpath(genpath(pwd))

load('ReachData_Bega.mat')
load('Q_Bega.mat')

%% global variable initialization

%mimimum volume to be considered for mobilization
global roundpar;
roundpar = 0;
   
%% Sediment classes definition
   
sed_range = [-5.5 , 3.5]; %range of sediment sizes considered in the model
class_size = 1.5; %amplitude of the sediment classes

global psi
psi =  sed_range(1):class_size:sed_range(2);   

clear sed_range class_size

%% Preprocessing
Network = graph_preprocessing(ReachData);

%% define external sediment inputs

clear Qbi_input; Qbi_input  = cell(timescale,1); [Qbi_input{:}] = deal(zeros(length(Network.II), length(psi)));
%[Qbi_input{2:end}] = deal([ones(1, length(psi) ).* 100.* eq50_Fi(1,:) ; zeros(length(NH)-1, length(psi))]);

%% define input sediment load in the deposit layer

deposit = [ReachData.deposit].*[ReachData.Length];
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

timescale = 200;

[data_output,extended_output] = DCASCADE_main( ReachData , Network , Q , timescale , Qbi_dep_in , Qbi_input  );


