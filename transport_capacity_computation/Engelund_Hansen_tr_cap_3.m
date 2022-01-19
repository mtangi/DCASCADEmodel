function [ Qtr_cap, pci ] = Engelund_Hansen_tr_cap_3(Fi_r_reach , Slope_reach , Wac_reach, v_reach, h_reach, varargin)

%ENGELUND_HANSEN_TR_CAP returns the value of the transport capacity (in Kg/s) 
%for each sediment class in the reach measured using the Engelund and Hansen equations

%% references
%Engelund, F., and E. Hansen (1967), A Monograph on Sediment Transport in Alluvial Streams, Tekniskforlag, Copenhagen.

%% read input

global psi

def_par =  'Molinas' ;

p = inputParser;
addOptional(p,'partitioning',def_par);
parse(p,varargin{:})

par_method = p.Results.partitioning ;

%% Transport capacity from Engelund-Hansen equations using the Molinas transport capacity fraction approach  (TCF, Molinas and Wu, 2000) 

if strcmpi(par_method , 'molinas')
 
    % find D values 
    D_values = [16 50 84];
    [D_changes] = D_finder_3(Fi_r_reach, D_values );

    D50 = D_changes(2);

    rho_s = 2650; % sediment densit [kg/m^3]
    rho_w = 1000; % water density [kg/m^3]
    g = 9.81;

    %friction factor
    C = (2*g.*Slope_reach.*h_reach)./(v_reach).^2;
    %dimensionless shear stress
    tauEH = (Slope_reach.*h_reach)./((rho_s/rho_w-1).*D50);
    %dimensionless transport capacity
    qEH = 0.05./C.* (tauEH).^(5/2);
    %dimensionful transport capacity m3/s 
    qEH_dim = qEH.*sqrt((rho_s./rho_w-1)*g*(D50).^3);%m3/s (%formula from the original cascade paper)

    QS_kg = qEH_dim.*Wac_reach.*rho_s; %kg/s

    QS_EH = QS_kg./rho_s; %m3/s

    %then the different sediment transport capacities have to be
    %splitted according to Molinas and saved into the Qbi_tr in
    %order to get the right structure for outputs.

    pci = Molinas_3( Fi_r_reach, h_reach, v_reach, Slope_reach, D_changes);
    Qtr_cap = pci.*QS_EH;
end

%% Transport capacity from Engelund-Hansen equations using bed material fraction approach (BMF, Molinas and Wu, 2000)

if strcmpi(par_method , 'bmf')
    
    dmi = 2.^(-psi)./1000; %sediment classes diameter (m)

    rho_s = 2650; % sediment densit [kg/m^3]
    rho_w = 1000; % water density [kg/m^3]
    g = 9.81;

    for d=1:length(dmi)
        
        %friction factor
        C = (2*g.*Slope_reach.*h_reach)./(v_reach).^2;
        %dimensionless shear stress
        tauEH = (Slope_reach.*h_reach)./((rho_s/rho_w-1).*dmi');
        %dimensionless transport capacity
        qEH = 0.05./C.* (tauEH).^(5/2);
        %dimensionful transport capacity m3/s 
        qEH_dim = qEH.*sqrt((rho_s./rho_w-1)*g*(dmi').^3);%m3/s (%formula from the original cascade paper)

        QS_kg = qEH_dim.*Wac_reach.*rho_s; %kg/s

        QS_EH = QS_kg./rho_s; %m3/s
        
        Qtr_cap = QS_EH.*Fi_r_reach;
    end
    
end

end
