function [roundvector] = rounddown(vector, roundpar)
%ROUNDDOWN rounds to the lowest decimal number given the approximation
%factor roundpar
% ex roundpar = 1  217.274 --->  213.2
%    roundpar = -1  217.274 --->  210

%%

if roundpar == 0
    roundvector = floor(vector);
elseif roundpar > 1 
    roundvector = floor(vector.*(10^roundpar))./(10^roundpar);
else
    roundvector = floor(vector./(10^(-roundpar))).*(10^(-roundpar));
end

end

