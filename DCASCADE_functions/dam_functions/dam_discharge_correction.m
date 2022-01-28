function [Q_new] = dam_discharge_correction(dam_considered, DamDatabase_active ,Network, release ,Q_t)
%dam_discharge_correction measures the change in discharge downstream
%each dam given the discharge of the dam and the locations of downstream
%dams

%% define dam order for the loop

%to determine the discharge for all reaches, i have to
%process the dams located upstream first
DamDatabase_considered = DamDatabase_active(dam_considered);
NH = Network.NH;
dam_order = NH(logical(sum(NH == [DamDatabase_considered.Node_ID]',1)));
outlet = Network.NH(end);

%% loop for each dam
Q_new = Q_t;

for d=1:length(DamDatabase_considered)
    
    %extract values for the dam
    dam_pos = find([DamDatabase_considered.Node_ID]==dam_order(d)); %position of the selected dam in DamDatabase
    
    Node_ID = [DamDatabase_considered(dam_pos).Node_ID];
    
    downstream_reaches = Network.Downstream.Path{Node_ID}{outlet}; %ID of the reaches downstream the dam
    
    node_ID_down = [DamDatabase_active.Node_ID]';
    node_ID_down(node_ID_down == Node_ID) = []; 
    [~,  down_pos] = intersect(downstream_reaches, node_ID_down); %position of the dams in the downstream path from d to outlet
    
    %derive discharge increase downstream of dam
    Q_diff = zeros(size(Q_t)); %Q_diff contains the difference of discharge between the dam reach and the downstram reaches. It forms a baseline discharge that will be added to the release
    Q_diff( downstream_reaches ) = max(Q_t(downstream_reaches) - Q_t(Node_ID),0);
        
    % find if there are flooded reaches belonging to other dams downstream d
    if isempty(down_pos)
        Q_new(downstream_reaches) = release(dam_pos) + Q_diff(:,downstream_reaches);
    else %if there are, the release is propagated downstream until the first barrier.
        Q_new(downstream_reaches(1:min(down_pos)) ) = release(dam_pos) + Q_diff(:,downstream_reaches(1:min(down_pos)));
    end
        
end


end
