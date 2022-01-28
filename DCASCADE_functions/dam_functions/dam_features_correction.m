function [  ReservoirHydraulics ] = dam_features_correction( DamDatabase_active , ReachData, Network , ResVolume_t , Slope_t , Node_el_t , Q_t )
%'dam_features_correction' changes the Slope and Width in flooded reaches

% input:
%      - 'DamDatabase' 1x1 struct containing the fields:
%           .id_FN : ID of the node where the dam is positioned
%           .FSL : elevation of the full supply level of the reservoir, sum of the elevation of id_FN and the dam heigth
%           .Name : cell array containing the name of the dams
%           .portfolio : vector containing '1' if the dam is present, '0' if absent
%      -'ReachData' Matrix with column variables defined below by ID_* 
%      -'Network' 1x1 struct containing for each node info on upstream and downstream nodes
%      -'calcResArea' vector containing for each dam the reservoir area [m2]

% output:
%      - 'ReservoirHydraulics' table containing the the features of the reservoirs on the river

%% 
% parameters for the "entire" reservoir
reservoirFromN=zeros(size(DamDatabase_active,1),1); % From-nodes of  reservoirs 
WL =zeros(size(DamDatabase_active,1),1); % Water level  of  reservoirs [m]
WE =zeros(size(DamDatabase_active,1),1); % Water elevation  of  reservoirs [masl]
Ares =zeros(size(DamDatabase_active,1),1); % Innundated Area [km2]
WRes_max=zeros(size(DamDatabase_active,1),1); % Reservoir width at dam site [km2]
maxLres=zeros(size(DamDatabase_active,1),1); % Reservoir length [km]

% parameters for specific nodes within the reservoir
ResInundNodes=cell(length(DamDatabase_active),1); % Innundated upstream nodes
W_reach=cell(length(DamDatabase_active),1); % Mean width in each innundated edge [m]
h_reach=cell(length(DamDatabase_active),1); % Mean water level in each innundated edge [m]
A_reach=cell(length(DamDatabase_active),1); % inudated surface at each innundated edge [m2]
Vl_reach=cell(length(DamDatabase_active),1); % Storage volume at each innundated edge [m2]
avg_flow_v=cell(length(DamDatabase_active),1); % average flow velocity at each innundated edge  [m/s]
Energy_slope_reach=cell(length(DamDatabase_active),1); % Active slope at each innundated edge 
c_res=zeros(length(DamDatabase_active),1); % Scaling factor reservoir width vs. Length

for d=1:length([DamDatabase_active.Node_ID])
    % definitions 
    
    reservoirFromN(d,1)=DamDatabase_active(d).Node_ID;
    % reservoirFromN(dam_i,1)=AggData(reservoirFromN(dam_i,1),ID_ToN);
    
    %find inverted WL-Volume function to find WL
    WL(d,1) = Reservoir_LSWConversion(ResVolume_t(d) ,DamDatabase_active(d).WL_table, 3, 1); %Water elevation level given the Reservoir volume;
    WE(d,1) =  WL(d,1) + Node_el_t( ReachData(reservoirFromN(d,1)).FromN) ; % Water elevation [masl]

    Ares(d,1) = Reservoir_LSWConversion(ResVolume_t(d) ,DamDatabase_active(d).WL_table, 3, 2); % if no area at FSL is given, load it from the model results
    %Ares(d,1) = DamDatabase(d).ResFSL_polyline(FSE(d,1)); %end % if no area at FSL is given, load it from the model results

    % find all nodes that are within the reservoir
    % all upstream nodes of the reservoir
    ResInundNodes{d,1}=find(cellfun(@length,Network.Upstream.Path{reservoirFromN(d)})>1); % all nodes upstream of a dam 

    % all nodes upstream of a dam below the average supply level, i.e., in the impoundment
    ResInundNodes{d,1}=ResInundNodes{d,1}( Node_el_t( [ReachData(ResInundNodes{d}).ToN] ) < WE(d) );

    % add this node to the reservoir nodes. This is important to get the
    % length of the reservoir right. 
    ResInundNodes{d,1}=unique(ResInundNodes{d,1}); 

    % if there are no upstream nodes: Add only the immediately upstream reach. 
    if isempty(ResInundNodes{d})
        ResInundNodes{d,1}=reservoirFromN(d);
    end 

    % calculate reservoir length:

    Lres=nan(size(ResInundNodes{d})); % distance of nodes from reservoir [m]

    for iii=1:length(ResInundNodes{d})
        Lres(iii)=Network.Downstream.Distance{ResInundNodes{d}(iii)}(reservoirFromN(d)); % distance of each innundated node from the dam[m]  
    end 

    if length(Lres)==1 % if there is only one innundated edges: assume that  the dam is located in the middle of the the current edge, and that half of the edge is innundated
        Lres(1)=ReachData(reservoirFromN(d,1)).Length;
    end  

    % check if some upstream edges are not fully inundated -> correct for
    % that effect. 

     A=WE(d)-Node_el_t( [ReachData(ResInundNodes{d}).ToN] ) ; % difference FSL-reservoir bottom upstream in an edge [m/m]
     B=WE(d)-Node_el_t( [ReachData(ResInundNodes{d}).FromN] ) ; % difference FSL-reservoir bottom downstream in an edge [m/m]

     B(B<0)=0; 

     delta_h=A-B; % difference in reservoir bottom elevation over an edge [m] 
     L=delta_h./Slope_t(ResInundNodes{d}); % innundated length (for edges that are fully innundated L equals the reach length). 

     if length(Lres)>1 % check if there is more than on innundated node
     L(ResInundNodes{d}==reservoirFromN(d))=0; % distance in node where reservoir is located is 0. 
     end 

     % if B == 0, a reach is not fully inundated: subtract not innundated length of that edge  
     if any(B==0)
     Lres(B<=0)=Lres(B<=0)-([ReachData(ResInundNodes{d}(B<=0)).Length]-L(B<=0)); 
     end 

    maxLres(d,1)=max(Lres); % MAximum path length in the reservoir  

    % calculate width at dam [m], and width regression     
    WRes_max(d,1)=2*Ares(d).*1E6/maxLres(d,1); % calculate reservoir width at the dam [ m ]. Assuming that the impoundment is of triangular shape
    c_res(d,1)=WRes_max(d)/(maxLres(d,1)); % decrease of reservoir width with increasing distance from the dam [m width/m length]

    % calculate mean width in each innundated reach.    

    % Mean distance in from dam in each reach [m]
    DistDam_inNodeii=Lres-0.5*L; % mean distance of each edge from the dam. 

    %  Width in each reach [m]: 
    W_reach{d}=WRes_max(d,1)-c_res(d,1)*DistDam_inNodeii; 

    % Area in each reach [m2]
    A_reach{d}=W_reach{d}.*L; 

    % mean water level in each reach [m] 
    h_reach{d}=mean([A;B],1);

    % stored volume in each reach 
    Vl_reach{d}=A_reach{d}.*h_reach{d}; 

    % average flow velocity and active slope in each reach 
    avg_flow_v{d} = DamDatabase_active(d).mean_inflow./(h_reach{d}.*W_reach{d});
    %avg_flow_v{d} = Q_t([ReachData(ResInundNodes{d}).reach_id])./(h_reach{d}.*W_reach{d});
    Energy_slope_reach{d} = ((avg_flow_v{d}).^2./(2*9.807));
    
end

%%%  write reservoir hydraulics table
ReservoirHydraulics = table(reservoirFromN,WE,WL,maxLres,WRes_max,c_res,ResInundNodes,Energy_slope_reach,W_reach,A_reach,h_reach,Vl_reach, avg_flow_v,'RowNames',{DamDatabase_active.Name}'); 
   
end
