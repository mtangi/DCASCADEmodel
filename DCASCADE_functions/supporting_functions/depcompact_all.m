function [Qbi_dep_comp] = depcompact_all(Qbi_dep, ReachData, Network , varargin )
%DEPCOMPACT_ALL uses function depcompact to compact all the layers in the
%cells in Qbi_dep

global roundpar

p = inputParser;
addOptional(p,'h_buffer',0);
addOptional(p,'h_layer',0.1);
addOptional(p,'CompSingle',0);

parse(p,varargin{:})

h_buffer = p.Results.h_buffer ;
h_layer = p.Results.h_layer ;
CompSingle = p.Results.CompSingle;

%% initialization 
Wac = [ReachData.Wac];
Lngt = [ReachData.Length];

Qbi_dep_comp  = cell(size(Qbi_dep,1),length(Network.II));

%% 
wb = waitbar(1/size(Qbi_dep,1) * size(Qbi_dep,2), ['timestep 1/' num2str(size(Qbi_dep,1))]); %open waitbar
timerout = 1000000;

for t = 1:size(Qbi_dep,1)

        
    time1   = clock;
    waitbar(t/size(Qbi_dep,1), wb, ['timestep ' num2str(t) '/' num2str(size(Qbi_dep,1)) '  -  ' num2str(ceil(timerout * (size(Qbi_dep,1) - t )) ) ' sec left' ]); % update waitbar


    for n = 1:  size(Qbi_dep,2)   
        
        Qbi_dep_comp{t,n} = depcompact(Qbi_dep{t,n}, Network.Upstream.distancelist{n} , Wac(n), Lngt(n) , h_buffer , h_layer);
        % start extraction
        if CompSingle ==1
            Qbi_dep_comp{t,n} = single(Qbi_dep_comp{t,n});
        end
        
        if any(round( sum(Qbi_dep_comp{t,n}(:,2:end),1), roundpar)  ~= round( sum(Qbi_dep{t,n}(:,2:end),1), roundpar))
            warning(['approximation error in reach ' num2str(n) ' at timestep ' num2str(t) ' = ' num2str(round( sum(Qbi_dep_comp{t,n}(:,2:end),1), roundpar)  - round( sum(Qbi_dep{t,n}(:,2:end),1), roundpar)) ])
        end

    end

    %measure time of routing
    time2   = clock;

    if mod(t , 10 ) == 0  %save time only at certain timesteps
         timerout = etime(time2, time1);
    end


end

close(wb); clear wb;

end