function out = split_cells(c_array)
%SPLIT_CELLS  Concatenate the contents of a cell array of row vectors.

for i= 1:numel(c_array)
    out{2*i - 1} = c_array{i}(1:round(length(c_array{i})./2));
    out{2*i} = c_array{i}(round(length(c_array{i})./2):end);
    
end