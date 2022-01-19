function [Qbi_sort] = sortdistance(Qbi, distancelist )
% SORTDISTACE sort the rows of the Qbi_incoming matrix by 
% increasing distance from the reach

%% code
[index,~] = find(  Qbi(:,1) == distancelist(distancelist ~= inf ) );

if ~isempty(index)
    
    Qbi_sort = Qbi(index,:);
    
else
    
    Qbi_sort = Qbi;
    
end
    
end
