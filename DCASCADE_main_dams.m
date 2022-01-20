function [data_output,extended_output,dam_output] = DCASCADE_main_dams( ReachData , Network , Q , timescale , DamDatabase, ORparameters, Qbi_dep_in , Qbi_input , dates_Q )
%% DCASCADE_main_dams
%
% INPUT :
%
% ReachData      = nx1 Struct defining the features of the network reaches
% Network        = 1x1 struct containing for each node info on upstream and downstream nodes
% Q              = txn matrix reporting the discharge for each timestep
% timescale      = length for the time horizion considered in the 
%
%
%----
% OUTPUT: 
%
% data_output      = struct collecting the main aggregated output matrices 
% extended_output  = struct collecting the raw D-CASCADE output datasets
% dam_outputs      = struct collecting main output for the reservoirs in the network

%% load global variables

global psi
global roundpar

%% Variables extraction from ReachData 

%reach properties
Lngt = [ReachData.Length];
n_Man = [ReachData.n];

Wac = repmat([ReachData.Wac],[ timescale 1 ]);

NH = Network.NH; %node hierarchy

outlet = Network.NH(end); % outlet reach ID identification
          
Q_original = Q;

%% variables definition 
    
%Variables to be calculated during the routing
clear Qbi_tr; Qbi_tr  = cell(timescale,1);     %Qbi_tr report the sediment mobilized present in the reach AFTER transfer
clear Qbi_mob; Qbi_mob  = cell(timescale,1);   %Qbi_mob report the sediment mobilized present in the reach BEFORE transfer

clear Qbi_dep; Qbi_dep  = cell(timescale,length(Network.II));

clear Fi_r_act; Fi_r_act = cell(timescale,1); %Fi_r_act contains the sed distribution of the active layer

clear Q_out ; Q_out = cell(timescale,1);

clear D50_AL; D50_AL = zeros(timescale,length(NH)); % support variable recording the D50 of the active layer in each reach in each timestep

%% parameters initialization

phi = 0.4; % sediment porosity in the active layer
minvel = 0.00001;

%% Slope and Node_el initialization 

clear Slope; Slope = zeros(timescale,length(NH)); 
clear Node_el ; Node_el = zeros(timescale,length(NH)+1);
clear EnergySlope; EnergySlope = zeros(timescale,length(NH)); 

% put a minimum value in the slope, to guarantee movement 
min_slope = min([ReachData.Slope]); %0.0005;

%initialize Slope
Slope(1,:) = max([ReachData.Slope],min_slope); 
Slope(2,:) = max([ReachData.Slope],min_slope); 

%initialize node elevation (for each reach the matrix reports the fromN elevation)
%the last column reports the outlet ToNode elevation (last node of the network), which can never change elevation.

Node_el(1,:) = [[ReachData.el_FN] [ReachData(outlet).el_TN]];
Node_el(2,:) = [[ReachData.el_FN] [ReachData(outlet).el_TN]];
Node_el(:,end) =  Node_el(1,end);
    
%% load Dam data
DamDatabase.portfolio; %attribute value of flow according to the percentile previously calculated
clear p

DamDatabase_active = DamDatabase([DamDatabase.portfolio] == 1) ;
ORparameters_active = ORparameters([DamDatabase.portfolio] == 1) ;

% find reaches in the reduced network with only flooded reaches and reaches
% downstream dams
[ReservoirHydraulics_full ] = dam_features_correction( DamDatabase_active , ReachData, Network, [DamDatabase_active.FSL_ResVolume], Slope(1,:) , Node_el(1,:) , Q(1,:));

flooded_reaches_full = ReservoirHydraulics_full.ResInundNodes;

%% initialize dam variables
clear start_vol
for d=1:length(DamDatabase_active) 
    start_vol(d) = Reservoir_LSWConversion(DamDatabase_active(d).WL_start ,DamDatabase_active(d).WL_table, 1, 3);
end

%reservoir stored volume [m3]
ResVolume = zeros(timescale+1, length(DamDatabase_active));
ResVolume(1:2,:) = repmat(start_vol,2,1);

%reservoir release [m3/s]
release = zeros(timescale+1, length(DamDatabase_active));

%reservoir full supply level [m]
FSL_ResVolume = zeros(timescale+1, length(DamDatabase_active));
FSL_ResVolume(1:2,:) = repmat([DamDatabase_active.FSL_ResVolume],2,1);

%reservoir sediment storage [m3]
sed_storage = zeros(timescale+1, length(DamDatabase_active));

%operating rule ID
id_OR = 3;

%% variables initialization
    
Qbi_tr{1} = single(zeros([size(Network.II),length(psi)]));  
Qbi_mob{1} = single(zeros([size(Network.II),length(psi)])); 

Q_out{1} = zeros( length(Network.II), length(psi));

% initialize sediment deposit in the reaches
for n = NH 
    %If no inputs are defined, initialize deposit layer with a single cascade with no volume and GSD equal to 0
    if isempty(Qbi_dep_in{n})
        Qbi_dep{1,n} = single([n, zeros(1,length(psi)) ]);
        Fi_r_act{1}(:,n) = zeros(1,length(psi)) ;
        D50_AL(1,n) = 0;
    else
        Qbi_dep{1,n} = single([ repmat(n,[size(Qbi_dep_in{n},1),1] ) Qbi_dep_in{n}]) ;
        Fi_r_act{1}(:,n) = sum(Qbi_dep_in{n},1)./sum(Qbi_dep_in{n},'all') ;
        D50_AL(1,n) = D_finder(Fi_r_act{1}(:,n));
    end
    
    Qbi_dep{2,n} = Qbi_dep{1,n};
end  
    
%set limit for erosion in 1 timestep, given by the parameter mlim, that
%is the maximum depth in meter that can be eroded
mlim = ones(1,length(ReachData))  * 0.1 ; 
V_lim_tot = round (mlim .* Wac .* Lngt , roundpar) ;
 
% support variable recording the total transport capacity in each reach in each timestep
tr_cap_sum = zeros(timescale, length(ReachData));
    
Qbi_tr{2} = Qbi_tr{1};

%% Routing scheme

tic   

%plot waitbar
wb = waitbar(1/timescale, ['timestep 1/' num2str(timescale)]); %open waitbar
timerout = 1000;

%start time loop
for t = 2: timescale-1
    
    %update waitbar
    time1   = clock;
    waitbar(t/timescale, wb, ['timestep ' num2str(t) '/' num2str(timescale) '  -  ' num2str(ceil(timerout * (timescale-t))/60) ' min left' ]); % update waitbar

        % change reservoir FSL according to reservoir sediment storage
    [FSL_ResVolume(t,:), sed_storage(t,:)] = changeFSLResVolume(DamDatabase_active,  Qbi_tr{t} , Qbi_dep(t,:), flooded_reaches_full,'phi',phi);

    %calculate dam mass balance and dam release
    [release(t,:), ResVolume(t+1,:), Q(t,:)] = dam_mass_balance(DamDatabase_active,Network, ResVolume(t,:), Q_original(t,:) ,'OperatingRule',id_OR,'ORparameters',ORparameters_active, 'date',dates_Q(:,t) , 'Q_input_estimates' , Q_original(t-1,[DamDatabase_active.Node_ID]) , 'FSL_ResVolume', FSL_ResVolume(t,:) );
    
    %calculate reach features for the flooded reaches
    [ReservoirHydraulics ] = dam_features_correction( DamDatabase_active , ReachData, Network, ResVolume(t,:) , Slope(t,:) , Node_el(t,:) , Q(t,:));
    
    % extract values of ReservoirHydraulics

    f_reach = cell2mat([ReservoirHydraulics.ResInundNodes]'); %inundate reaches ID
    EnergySlope(t,:) = Slope(t,:); %in all reaches but the flooded ones, energy slope is equal to the slope 
    EnergySlope(t,f_reach) = min( cell2mat([ReservoirHydraulics.Energy_slope_reach]'), Slope(t,f_reach) );
    Wac(t,f_reach) = max( cell2mat([ReservoirHydraulics.W_reach]'), [ReachData(f_reach).Wac] );
    
    %calculate new water dept for all reaches
    %via Manning equation
    h = (Q(t,:).*n_Man./(Wac(t,:).*sqrt( EnergySlope(t,:) ))).^(3/5);
    h(f_reach) = max( cell2mat([ReservoirHydraulics.h_reach]'), h(f_reach)); % water depth for flooded reaches

    v = 1./n_Man.*h.^(2/3).*sqrt( EnergySlope(t,:) ); 
    v(f_reach) = min( cell2mat([ReservoirHydraulics.avg_flow_v]') , v(f_reach)); % water velocity for flooded reaches
      
    %variables initialization for the timestep
    Qbi_tr{t+1} = single(zeros( size(Qbi_tr{1}) ));
    Qbi_mob{t} = single(zeros( size(Qbi_mob{1}) ));
    Q_out{t} = zeros( size(Q_out{1}) );
            
    %loop for all reaches
    for n = NH

        %reach_specific state variables initialization
        Qbi_tr{t+1}(:,n,:) = zeros( size(Qbi_tr{t}(:,n,:)) );
        Qbi_dep{t+1,n} = single([n, zeros(1,length(psi)) ]);
        
        %output variable initialization
        Qbi_mob{t}(:,n,:) = zeros( size(Qbi_mob{t-1}(:,n,:)) );
        Q_out{t}(n,:) = zeros( size( Q_out{t-1}(n,:)) );

        V_dep_old = double(Qbi_dep{t,n}); % extract the deposit layer of the reach from the relative cell in the previous timestep

        %%% 1) extract the deposit layer from the storage matri and load the incoming cascades

        Qbi_incoming = [(1:length(NH))' squeeze(Qbi_tr{t}(:, n,:)) ; n Qbi_input{t}(n,:) ]; 
        Qbi_incoming(sum(Qbi_incoming(:,2:end),2)==0,:) = [];

        if isempty(Qbi_incoming); Qbi_incoming = [n zeros(1,length(psi))]; end %put an empty cascade if no incoming volumes are present (for computation)

        %sort incoming matrix accoring to distance, in this way
        %sediment coming from closer reaches will be deposited first 
        [Qbi_incoming] = sortdistance(Qbi_incoming, Network.Upstream.distancelist{n} );

        %%% 2) find cascades to be included into the active layer according to the limit V_lim_tot, and use the cumulative GSD to compute tr_cap

        %find the volume of sediment from the incoming load (V_inc2act) and deposit layer (V_dep2act) to be included in the active layer
        [V_inc2act , V_dep2act ,  V_dep , Fi_r_act{t}(:,n)] = layer_search (Qbi_incoming, V_dep_old, V_lim_tot(n));

        %if i have almost no sediment in my reach (less then 1/100 of the active layer limit, i put the sed distribution to the values in the previous timestep)
%             if sum(Qbi_incoming(:,2:end),'all') + sum(V_dep_old(:,2:end),'all') < V_lim_tot(n)/100
%                 Fi_r_act{t}(:,n) = Fi_r_act{t-1}(:,n);
%             end

        if sum(Fi_r_act{t}(:,n))==0; Fi_r_act{t}(:,n) = Fi_r_act{t-1}(:,n); end % in case the active layer is empty, i use the GSD of the previous timesteep

        D50_AL(t,n) = D_finder(Fi_r_act{t}(:,n));

        %calculate transport capacity using the Fi of the active layer, the resulting tr_cap is in m3/s and is converted in m3/day
        tr_cap = Engelund_Hansen_tr_cap_3(Fi_r_act{t}(:,n) , EnergySlope(t,n)  , Wac(t,n), v(n) , h(n), 'partitioning','bmf' )' .* 24.*60.*60;

        tr_cap_sum(t,n) = sum(tr_cap);
        
        %%% 3) Deposit the cascades in the active layer until the volume mobilized for each class is equal to the tr_cap

        if sum(tr_cap) < ( sum(V_dep2act(:,2:end),'all') + sum(V_inc2act(:,2:end),'all') ) %if the total transport capacity is lower than the active layer volume...
              %... deposit part of the active layer cascades, 
              %    proportionally to their volume and the volume of the active layer

              [V_mob, V_dep ] = tr_cap_deposit( V_inc2act, V_dep2act, V_dep, tr_cap);      

        else
            % if not, the mobilized layer is equal to the active
            % layer
            V_mob = matrix_compact([V_dep2act ; V_inc2act]);

        end 

        % (after this passage, V_mob contains only the volumes actually mobilized)     
        Qbi_dep{t+1,n} = single(V_dep);

        %remove empty rows
        Qbi_dep{t+1,n}(sum(Qbi_dep{t+1,n}(:,2:end),2)==0,:)=[];

        % Qbi_mob contains the volume mobilized in the reach, that is about to be transfer downstream
        Qbi_mob{t}(V_mob(:,1),n,:) = single(V_mob(:,2:end)) ; 

        % if removing empty rows leaves only an Qbi_dep{t,n} empty
        % matrix, put an empty layer
        if isempty( Qbi_dep{t+1,n}); Qbi_dep{t+1,n} = single([n zeros(1,length(psi))]); end

        %%% 4) Compute the changes in bed elevation

        % modify bed elevation according to increased deposit

        Delta_V = sum(Qbi_dep{t+1,n}(:,2:end),'all') -  sum(Qbi_dep{t,n}(:,2:end),'all');
        Node_el(t+1,n) = Node_el(t,n) + Delta_V/( ( sum(Wac(t,[n Network.Upstream.Node{n}]) .* Lngt([n Network.Upstream.Node{n}]) ) ) * (1-phi) );

    %(end of the reach loop)
    end

    %%% 5) Move the mobilized volumes to the destination reaches according to the sediment velocity

    %loop for all reaches, now that i have the Fi_r and thus can compute transfer rates for all reaches
    clear Qbi_tr_t
    for n = NH
        
        %load mobilized volume for reach n
        V_mob = zeros(length(NH),length(psi)+1);
        V_mob(:,1) = 1:length(NH);

        V_mob(:,2:length(psi)+1) = squeeze(Qbi_mob{t}(:,n,:));
        V_mob = matrix_compact(V_mob);
        
        %calculate GSD of mobilized volume
        Fi_mob = sum(V_mob(:,2:end),1)'./sum(V_mob(:,2:end),'all');
        if isnan(Fi_mob); Fi_mob = Fi_r_act{t}(:,n);           end

        %calculate sediment velocity for the mobilized volume in each reach
        [ v_sed ] = velocity_EH( repmat( Fi_mob,[1 length(NH)] ) , EnergySlope(t,:) , Wac(t,:) , v , h , 'minvel', minvel,'phi',phi,'IDformula', 1 ) ;

        %transfer the sediment volume downstream
        [Qbi_tr_t, Q_out_t  ] = sed_transfer_simple( V_mob , n , v_sed.*(60*60*24) , Lngt, Network );

        % Sum the volumes transported from reach n with all the other
        % volumes mobilized by all the other reaches at time t

        Qbi_tr{t+1} = Qbi_tr{t+1} + single(Qbi_tr_t);
        Q_out{t} =  Q_out{t} + Q_out_t;

    end

    % change the slope accordingly to the bed elevation
    [Slope(t+1,:), Node_el(t+1,:)] = change_slope( Node_el(t+1,:) ,Lngt, Network, min_slope );
            
    %measure time of routing
    time2   = clock;

    if (mod(10, t) == 0)   %save time only at certain timesteps
         timerout = etime(time2, time1);
    end
    
%end of the time loop  
end

close(wb); clear wb time1 time2 timerout Qbi_tr_t Q_out_t f_reach
toc 

%% output processing

% aggregated matrixes
QB_mob = cellfun( @(x)double(sum(x,3)),Qbi_mob(1:timescale-1),'UniformOutput',0); %total sediment mobilized in each reach (column), divided by reach provenance (row)
QB_mob_sum = cell2mat( cellfun(@(x)double(sum(x,[1,3])),Qbi_mob(1:timescale-1),'UniformOutput',0)); %total sediment mobilized in each reach (column)

%total sediment delivered in each reach (column), divided by reach provenance (row)
QB_tr = cellfun( @(x)double(sum(x,3)),Qbi_tr(1:timescale-1),'UniformOutput',0); 

%total material in the deposit layer
V_dep = cell2mat(cellfun( @(x)double(sum(x(:,2:end),'all')) , Qbi_dep(1:timescale-1,:),'UniformOutput',0));

% total volume in the deposit layer for each timestep, divided by sed.class
V_class_dep = cellfun( @(x)double(sum(x(:,2:end),1)) ,Qbi_dep(1:timescale-1,:),'UniformOutput',0); 

%total material in a reach in each timestep (both in the deposit layer and mobilized layer)
tot_sed = cell2mat(cellfun( @(x)double(sum(x(:,2:end),'all')) , Qbi_dep(1:timescale-1,:),'UniformOutput',0)) + cell2mat(cellfun(@(x)sum(x,1), QB_tr(1:timescale-1), 'UniformOutput', false)) ;

%sediment mobilized in each reach divided by class 
Qbi_mob_class = cellfun(@(x)squeeze(double(sum(x,1))),Qbi_mob(1:timescale-1),'UniformOutput',0);

% D50 of the material in a reach for each timestep 
D50_tot = zeros(size(Qbi_mob_class,1), size(V_class_dep,2));

%total material in a reach in each timestep, divided by class 
tot_sed_class = Qbi_mob_class;

for t = 1:length(Qbi_mob_class)
    tot_sed_class{t} = Qbi_mob_class{t} + reshape(cell2mat(V_class_dep(t,:)),[length(psi) , size(V_class_dep,2) ])';
    Fi_tot_t = tot_sed_class{t} ./ sum(tot_sed_class{t},2);
    Fi_tot_t(isnan(Fi_tot_t)) = 0;
    for i=1: size(V_class_dep,2)
        D50_tot(t,i) = D_finder(Fi_tot_t(i,:)');
    end
end

%total material in a reach in each timestep, divided by class 
tot_sed_class = cell(size(psi));

for c = 1:length(tot_sed_class)
    tot_sed_class{c} = cell2mat(cellfun( @(x)double(sum(x(:,c+1),'all')) , Qbi_dep(1:timescale-1,:),'UniformOutput',0)) +  cell2mat( cellfun( @(x)sum(x(:,:,c), 1) , Qbi_tr(1:timescale-1) , 'UniformOutput', false) );
end

%total sediment volume leaving the network
outcum_tot = cell2mat(cellfun(@(x) sum(x,'all'), Q_out(1:timescale-1), 'UniformOutput',0));

% set all NaN transport capacity to 0;
tr_cap_sum(isnan(tr_cap_sum)) = 0;

% set all NaN active layer D50 to 0;
D50_AL(isnan(D50_AL)) = 0;

%% compress the Qbi_dep matrix

%[Qbi_dep] = depcompact_all(Qbi_dep, ReachData, Network,'h_layer',1 , 'CompSingle',1);

   %% output struct definition
% data_plot contais the most important D_CASCADE outputs

data_output = cell(1,2);

data_output{1,1} = 'Wac';
data_output{2,1} = 'Mobilized volume';
data_output{3,1} = 'Total sed in the reach';
data_output{4,1} = 'Slope';
data_output{5,1} = 'D50 tot';
data_output{6,1} = 'D50 active layer';
data_output{7,1} = 'Total trasport capacity';
data_output{8,1} = 'Deposited volume';
data_output{9,1} = 'Discharge';
data_output{10,1} = 'Total sed in the reach - per class';

data_output{1,2} = repmat(Wac,[size(Qbi_dep,1),1]); 
data_output{2,2} = QB_mob_sum;
data_output{3,2} = tot_sed;
data_output{4,2} = Slope;
data_output{5,2} = D50_tot; 
data_output{6,2} = D50_AL; 
data_output{7,2} = tr_cap_sum;
data_output{8,2} = V_dep;
data_output{9,2} = Q(1:timescale,:);
data_output{10,2} = tot_sed_class;

% all other outputs are included in the extended_output cell variable
extended_output = cell(1,2);

extended_output{1,2} = Qbi_tr; 
extended_output{2,2} = Qbi_mob;
extended_output{3,2} = Q_out;
extended_output{4,2} = Qbi_dep;
extended_output{5,2} = Fi_r_act; 
extended_output{6,2} = D50_AL;
extended_output{7,2} = Node_el;

extended_output{1,1} = 'Qbi_tr'; 
extended_output{2,1} = 'Qbi_mob';
extended_output{3,1} = 'Q_out';
extended_output{4,1} = 'Qbi_dep';
extended_output{5,1} = 'Fi_r_ac';
extended_output{6,1} = 'D50_AL';
extended_output{7,1} = 'Fi_yield';
extended_output{8,1} = 'Node_el';

%%
%all the dams output are contained in dam_output

dam_output = cell(1,2);

dam_output{1,2} = DamDatabase_active; 
dam_output{2,2} = ORparameters_active; 
dam_output{3,2} = ResVolume;
dam_output{4,2} = release;
dam_output{5,2} = FSL_ResVolume;
dam_output{6,2} = sed_storage; 
dam_output{7,2} = EnergySlope;
dam_output{8,2} = flooded_reaches_full;

dam_output{1,1} = 'DamDatabase_active'; 
dam_output{2,1} = 'ORparameters_active';
dam_output{3,1} = 'ResVolume';
dam_output{4,1} = 'release';
dam_output{5,1} = 'FSL_ResVolume';
dam_output{6,1} = 'sed_storage';
dam_output{7,1} = 'EnergySlope';
dam_output{8,1} = 'flooded_reaches_full';

end