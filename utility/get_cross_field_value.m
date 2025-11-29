function subS = get_cross_field_value(S,i)
% GET_CROSS_FIELD_VALUE
%
%   Summary: Extracts the element at index 'i' from every field in the input
%   scalar structure 'S', and creates a new scalar structure 'subS' containing
%   only those extracted elements.
%
%   Inputs:
%       S (struct): A scalar structure whose fields contain arrays (vectors,
%                   cell arrays, etc.).
%       i (numeric): The single index to extract from each field array.
%
%   Output:
%       subS (struct): A new scalar structure containing the i-th element
%                      from each corresponding field in S.
%
%   Copyright (c) 2014-2025 Christian Paolo Richter
%   University of Osnabrueck

% Get a cell array 'fn' containing the names of all fields in the input structure S.
fn = fieldnames(S);

% Begin loop: Iterate through each field name obtained in 'fn'.
% 'idxField' tracks the current field's position in the list 'fn'.
for idxField = 1:numel(fn)

    % Dynamic Field Access and Extraction:
    % 1. (fn{idxField}): Gets the field name string (e.g., 'Name').
    % 2. S.(fn{idxField}): Accesses the array content of that field in the
    %    input structure S (e.g., S.Name, which might be a vector).
    % 3. S.(fn{idxField})(i): Indexes into that array content and extracts
    %    the element at position 'i'.
    % 4. subS.(fn{idxField}) = ...: Assigns the single extracted element
    %    to a new field with the same name in the output structure 'subS'.
    subS.(fn{idxField}) = S.(fn{idxField})(i);

end %for
end %fun