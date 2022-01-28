function [release_t, ResVolume_post, Q_t ] = dam_mass_balance(DamDatabase_active, Network,  ResVolume_pre, Q_t, varargin)
%DAM_MASS_BALANCE defines the dam reservoir volume and the release given
%the operating rule of the dam based on the inflow and volume at t. 

%% define dam order for the loop

%to determine the discharge for all reaches, i have to
%process the dams located upstream first

NH = Network.NH;
dam_order = NH(logical(sum(NH == [DamDatabase_active.Node_ID]',1)));

%% load additional input 

def_OR = ones(length(DamDatabase_active));
def_ORparam = [];

p = inputParser;
addOptional(p,'OperatingRule',def_OR);
addOptional(p,'ORparameters',def_ORparam);
addOptional(p,'date',0);
addOptional(p,'Q_input_estimates',0);
addOptional(p,'FSL_ResVolume',[DamDatabase_active.FSL_ResVolume]);

parse(p,varargin{:})

id_OR = p.Results.OperatingRule ;
ORparameters_active = p.Results.ORparameters ;
date_timestep = p.Results.date; %date relative to timestep t

% estimated discharge entering the reservoir at time t, used for release
% decision-making; (default: estimated inflow = inflow in the previous timestep)
Q_input_estimates = p.Results.Q_input_estimates;

FSL_ResVolume_t = p.Results.FSL_ResVolume; %reservoir volume at full supply level at time t

%% define operating rule

release_t = zeros(size(ResVolume_pre));
ResVolume_post = zeros(size(ResVolume_pre));

%for each dam
for d=1:length(DamDatabase_active)
    
    %extract values for the dam
    dam_pos = find([DamDatabase_active.Node_ID]==dam_order(d)); %position of the selected dam in DamDatabase

    % input discharge at time t
    Q_input = Q_t(DamDatabase_active(dam_pos).Node_ID);
            
    %if the reservoir volume exceed the FLS after the input, i define a
    %minimum release to avoid activating the spillways. Else, the minimum
    %release is defined by the dam feature min_discharge
    if ResVolume_pre(dam_pos) + Q_input *(24*60*60)/1E6 > FSL_ResVolume_t(dam_pos)
        minrelease =  max((ResVolume_pre(dam_pos) - FSL_ResVolume_t(dam_pos))*1E6/(24*60*60) + Q_input,DamDatabase_active(dam_pos).min_discharge);
    else
        minrelease = DamDatabase_active(dam_pos).min_discharge;
    end

    maxrelease = min( ResVolume_pre(dam_pos)*1E6/(24*60*60) + Q_input , DamDatabase_active(dam_pos).design_discharge);

    %if the reservoir volume exceed the FLS after the input and the maximum
    %possible release, the spillways are always activated and the release
    %is obtained as the difference between the incoming discharge and the
    %remaining volume left in the reservoir
    if ResVolume_pre(dam_pos) + ( Q_input - DamDatabase_active(dam_pos).design_discharge) *(24*60*60)/1E6 > FSL_ResVolume_t(dam_pos)
        release_t(dam_pos) = Q_input - ( FSL_ResVolume_t(dam_pos) - ResVolume_pre(dam_pos))/(24*60*60)*1E6;

    % if the reservoir is empty, I release the minimum discharge possible given the inflow in the timestep
    elseif ResVolume_pre(dam_pos) + ( Q_input - DamDatabase_active(dam_pos).min_discharge)*(24*60*60)/1E6  < 0 
        release_t(dam_pos) =  Q_input;   
        
    else

        switch id_OR

            % CASE 1: release equal to input
            % The dam releases exactly the input discharge, unless it is higher
            % then the design discharge or lower then the minimum discharge
            case 1
                release_t(dam_pos) = min( maxrelease, max(minrelease,Q_input ));

            % CASE 2: constant target
            % we try to keep the reservoir water level at a constant
            % value. Thus, if the reservoir level is in target, release is
            % equal to input, unless it exceeds the design discharge. The
            % input needed is the SL target (SL_target)
            case 2

                WL_target = ORparameters_active(dam_pos).WL_constant; %Surface level target
                ResVol_target =  Reservoir_LSWConversion(WL_target ,DamDatabase_active(dam_pos).WL_table, 1, 3) ; %Reservoir volume associated to the target
                Q_predict = Q_input_estimates(dam_pos); % prediction of input discharge, needed to define the release 
                opt_release = max(0 , ResVolume_pre(dam_pos) + Q_predict*(24*60*60)/1E6 - ResVol_target)*1E6/(24*60*60); %optimal release to meet the target
                release_t(dam_pos) = min( maxrelease , max(minrelease,opt_release ));

            % CASE 3: semestral target
            % the target reservoir level increases or decreases linearly
            % between two user-specified target and dates. Thus, the
            % reservoir level target changes daily. The needed input are
            % the current date for the timestep (date), and the dates and
            % values of the target maximum and minimum SL (SLtarget_MinMax)
            case 3

                WL_target = findWL_twotargets(ORparameters_active(dam_pos), date_timestep); % find target water level
                ResVol_target = Reservoir_LSWConversion(WL_target ,DamDatabase_active(dam_pos).WL_table, 1, 3);  %Reservoir volume associated to the target WL
                Q_predict = Q_input_estimates(dam_pos); % prediction of input discharge, needed to define the release;
                opt_release = max(0 , ResVolume_pre(dam_pos) + Q_predict*(24*60*60)/1E6 - ResVol_target)*1E6/(24*60*60); %optimal release to meet the target
                release_t(dam_pos) = min( maxrelease , max(minrelease,opt_release )); %actual release given the upper and lower release limits

            case 4

        end %end of the operating strategy CASE

    end %end of the spillway if

    %calculate hydrological changes downstream dams
    [Q_t] = dam_discharge_correction(dam_pos, DamDatabase_active , Network, release_t(dam_pos) , Q_t );

    %Volume at the end of the daily timestep is the difference
    %between incoming and released volume;
    ResVolume_post(dam_pos) = max(ResVolume_pre(dam_pos) + (Q_input - release_t(dam_pos))*(24*60*60)/1E6,0);
    
end

end

