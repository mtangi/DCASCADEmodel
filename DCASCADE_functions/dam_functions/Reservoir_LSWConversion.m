function [output] = Reservoir_LSWConversion(input ,input_table, id_input, id_output)
%Reservoir_LSWConversion interpolates the current value of water level(L),
%surface(S) or volume (V) of a reservoir (according to id_output), given some
%standard values of these parameters contained in table and an input value
%of a different parameter (defined by id_output).

% id_input and output
% 1: Supply level
% 2: Reservoir area
% 3: Reservoir volume

%% extract table values 
output = zeros(size(input));

for i=1:length(input)
table_input = input_table{:,id_input};
table_output = input_table{:,id_output};

pos_high = find(input(i)<table_input,1);
pos_low = pos_high - 1;

if sum(input(i)>table_input) >= length(table_input)
    warning('the input value is higher then the tabulated values');
end

%extract the output value by performing a linear interpolation between
%the lower and upper value of the input value, given by the table.

output(i) = table_output(pos_low)+( (input(i) - table_input(pos_low)) / (table_input(pos_high) - table_input(pos_low)) * (table_output(pos_high) - table_output(pos_low)) );

end

end

