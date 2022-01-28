function [ tr_cap ] = Wong_Parker_formula( D50, Slope, Wac, h)

%WONG_PARKER_TR_CAP returns the value of the transport capacity (in m3/s)
%for each sediment class in the reach measured using the Wong-Parker
%equations 

% This function is for use in the D-CASCADE toolbox

%% references
%Wong, M., and G. Parker (2006), Reanalysis and correction of bed-load relation of Meyer-Peter and M�uller using their own database, J. Hydraul. Eng., 132(11), 1159�1168, doi:10.1061/(ASCE)0733-9429(2006)132:11(1159).

%% Transport capacity from Wong-Parker equations
global psi

dmi = 2.^(-psi)./1000; %sediment classes diameter (m)
rho_s = 2600; % sediment densit [kg/m^3]
rho_w = 1000; % water density [kg/m^3]
g = 9.81;

%Wong_Parker parameters 

alpha = 3.97;
beta = 1.5;
tauC = 0.0495;

% alpha = 4.93;
% beta = 1.6;
% tauC = 0.0470;

%dimensionless shear stress
tauWP = (Slope*h)/((rho_s/rho_w-1)*D50);
%dimensionless transport capacity
qWP = alpha* (max(tauWP - tauC,0) )^(beta);
%dimensionful transport capacity m3/(s*m) 
qWP_dim = qWP * sqrt((rho_s/rho_w-1)* g * (D50)^3); %m3/(s*m) (formula from the original cascade paper)

QS_WP = qWP_dim * Wac; %m3/s

tr_cap = QS_WP; %m3/s

end