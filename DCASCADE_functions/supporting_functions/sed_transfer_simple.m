function [Qbi_tr_t, Q_out_t , setplace, setout] = sed_transfer_simple(V_mob , n , v_sed_day , Lngt, Network )
% SED_TRANSFER_SIMPLE takes the matrix of the mobilized layers(V_mob) and the vector of
% the sed velocity for each class(v_sed_id) in a reach (n) and returns the 3D matrices containing the
% sed volumes in V_mob in their new position at the end of the timestep.
%
% This simple version of the function represents the volume as a point
% sediment parcel delivered from the ToN of the reach n. Thus the volume
% has a single destination reach and it never get split.

%% initailize parameters
outlet = Network.NH(end);

%% find start and end reach of the sed volume after the timestep
    
    %reach_dest is the id of the reach where the sed. volume stops after the timestep 
    %p_dest is the position from the from_node of the id reach where the sed. volume stops after the timestep 

    if n == outlet  
       reach_dest = repmat( n , size(v_sed_day,1) ,1);
       p_dest = Lngt(n) + v_sed_day(:,n);    
    else   
        %to find p_end, i track the position of a sediment parcel starting
        %from the To_Node of the reach n (i.e. the From_node of the downstream reach).
        [p_dest , reach_dest] = track_sed_position( cell2mat(Network.Downstream.Node(n)) , v_sed_day , Lngt , Network ); 
        p_dest = p_dest + Lngt(n);
    end

    %downdist contains the distanche from the starting reach to all reaches
    %downstream
    downdist = Network.Downstream.Distance{n,1} ;
   
    %% find position of the sediment volume
        
    %setout is equal to 1 if the volume in the sed.class left the network
    %via the outlet
    setout = p_dest - Lngt(reach_dest)' - downdist(reach_dest)' > 0;
     
    %in each row, setout is equal to 1 in the reach where the sed. volume
    %of each class is delivered 
    setplace = zeros([size(v_sed_day,1) length(downdist) ]);
    setplace(sub2ind(size(setplace),1:size(v_sed_day,1),reach_dest') ) = 1 ;
    
    setplace(setout==1,:) = 0;

    %% place volume to destination reach 
    
    Qbi_tr_t = zeros (length(Lngt), length(Lngt) , size(setplace,1));
    Q_out_t = zeros (length(Lngt), size(setplace,1));
    
    for c = 1:size(setplace,1)
       
        Qbi_tr_t(V_mob(:,1),:,c) = V_mob(:,c+1) .* setplace(c,:);
        Q_out_t(V_mob(:,1),:) = V_mob(:,2:end) .* setout';
        
    end

    
end

