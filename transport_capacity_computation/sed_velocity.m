function [ v_sed ] = sed_velocity( indx_tr_cap , indx_partition,  Fi_r ,  Slope_reach , Q, Wac, v , h , varargin )
%VELOCITY_AW returns the velocity of the sediment (in m/s) for each sediment
%class for each reach using the Engelund Hansen equations (1967)
%
%OUTPUTS:
% 
% v_sed: [cxn] matrix reporting the velocity for each sediment class c for
%        each reach n [m/s]

%% optional input selection

def_minvel = 0;
def_phi = 0.4;
def_indx_velocity = 1;

p = inputParser;
addOptional(p,'minvel',def_minvel);
addOptional(p,'phi',def_phi);
addOptional(p,'velocity_formula',def_indx_velocity);

parse(p,varargin{:})

minvel = p.Results.minvel ;
phi = p.Results.phi ;
indx_velocity = p.Results.velocity_formula;

%% active layer definition

%active layer as 10% of the water column depth 
L_a = 0.1.*h ; %characteristic vertical length scale for transport.

%alternative: active layer as 2*D90 (Parker, 2008)
%L_a = 2*D_finder_3(Fi_r_reach, 90 );

%% D50 definition
global psi
dmi = 2.^(-psi)'./1000; %grain size classes[m]

% find D values 
D50 = D_finder(Fi_r, 50 );

%% sediment velocity with fractional trasport capacity

%  by measuring the trasport capacity for the single sed.classes
%  indipendenty, we obtain different values of sed. velocity

if indx_velocity == 1

    %choose transport capacity formula
    switch indx_tr_cap
        case 1
            tr_cap_formula = @(D50,Fi_r)Parker_Klingeman_formula( Fi_r, D50, Slope_reach, Wac , h);
        case 2
            tr_cap_formula = @(D50,Fi_r)Wilcock_Crowe_formula(Fi_r, D50, Slope_reach, Wac , h);
        case 3
            tr_cap_formula = @(D50)Engelund_Hansen_formula( D50 , Slope_reach , Wac, v , h );
        case 4
            tr_cap_formula = @(D50)Yang_formula( Fi_r, D50 , Slope_reach , Q, v, h );
        case 5
            tr_cap_formula = @(D50)Wong_Parker_formula( D50 ,Slope_reach, Wac ,h );
        case 6
            tr_cap_formula = @(D50)Ackers_White_formula( D50,  Slope_reach , Q, v, h);
    end

    if any(indx_tr_cap == [1,2]) % if I use fractional transport formulas...
        
        tr_cap = zeros(length(dmi), length(Slope_reach));
        % ... run the tr.cap function indipendently for each class, setting
        % the frequency of each class = 1
        for d=1:length(dmi)
            Fi_r = zeros(size(dmi));
            Fi_r(d) = 1;
            Fi_r = repmat(Fi_r,1,length(Slope_reach));
            
            tr_cap_class = tr_cap_formula(dmi(d), Fi_r);
            tr_cap(d,:) = tr_cap_class(d,:);
        end
        
    else
        tr_cap = cell2mat(arrayfun(tr_cap_formula,dmi,'UniformOutput',false));
    end
    
    %calculate velocity
    v_sed = max( tr_cap./( Wac .* L_a .* (1-phi) ) , minvel);
    
    v_sed(:,L_a==0) = minvel;
    
end

%% sediment velocity with total trasport capacity

% sediment velocity found in this way is constant for all sed.classes
if indx_velocity == 2
    [ Qtr_cap, pci ] = tr_cap_junction( indx_tr_cap , indx_partition , Fi_r , D50 ,  Slope_reach, Q, Wac, v , h );
    v_sed = max( Qtr_cap./( Wac .* L_a.*(1-phi) .* pci ) , minvel);
end


end