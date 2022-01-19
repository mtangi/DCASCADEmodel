function [Slope_t, Node_el_t] = change_slope(Node_el_t, Lngt, Network , min_slope)
%CHANGE_SLOPE modify the Slope vector according to the changing elevation of
%the nodes: It also guarantees that the slope is not negative or lower then
%the min_slope value by changing the node elevation bofore findin the SLlpe

%% define minimum reach slope
if nargin < 4
    min_slope = 0;
end

%% initialization

outlet = Network.NH(end);

down_node = Network.Downstream.Node;
down_node{outlet} = length(Node_el_t);
down_node = cell2mat(down_node);

Slope_t = zeros(size(Lngt));

%% loop for all reaches

for n=1:length(Lngt)
    
    % find the minimum node elevation to guarantee Slope > min_slope
    min_node_el = min_slope *  Lngt(n) + Node_el_t(down_node(n)); 

    % change the noide elevation if lower to min_node_el
    Node_el_t(n) = max(min_node_el, Node_el_t(n) );
    
    % find the new slope
    Slope_t(n) = (Node_el_t(n) - Node_el_t(down_node(n))) / Lngt(n);
    
end

end

