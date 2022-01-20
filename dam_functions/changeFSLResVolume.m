function [FSL_ResVolumenew , sed_storage] = changeFSLResVolume(DamDatabase_active,  Qbi_tr_t , Qbi_dep_t , flooded_reaches_full, varargin)
%changeFSLResVolume change the reservoir volume at full supply level according to the sediment
%storage in the reservoir;
def_phi = 0;

p = inputParser;
addOptional(p,'phi',def_phi);
parse(p,varargin{:})

phi = p.Results.phi ;

%%
FSL_ResVolumestart = [DamDatabase_active.FSL_ResVolume]; %reservoir volume at t [10^6 m3]
FSL_ResVolumenew = zeros(size(FSL_ResVolumestart)); %reservoir volume at t=1 [10^6 m3]
sed_storage = zeros(size(FSL_ResVolumestart)); %reservoir storage at t [m3]

for d=1:length(FSL_ResVolumestart)
    
    sed_storage_reach = cell2mat(cellfun(@(x)sum(x(:,2:end),'all'),Qbi_dep_t(flooded_reaches_full{d}),'UniformOutput',0))+sum(Qbi_tr_t(:,flooded_reaches_full{d},:),[1 3]);

    sed_storage(d) = sum(sed_storage_reach); % sediment storage in m3
    
    FSL_ResVolumenew(d) = FSL_ResVolumestart(d) - sed_storage(d)/(1-phi)/1E6;
        

end

end
