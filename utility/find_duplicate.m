function [V,I] = find_duplicate(x)
% FIND_DUPLICATE identifies duplicate values in a vector and returns the values
% and a logical matrix indicating their positions.
%
%   INPUTS:
%   x: vector (N x 1 or 1 x N); The input array (e.g., a list of cluster IDs).
%
%   OUTPUTS:
%   V: vector; The unique values found in x that occur more than once (the duplicates).
%   I: logical matrix (N x M); A matrix where N is the length of x, and M is the length of V.
%      I(i, j) is TRUE if x(i) equals V(j).
%
%   Copyright (c) 2018-2025 Christian Paolo Richter
%   University of Osnabrueck
%

UUID = unique(x);
% A is a logical matrix comparing every element of x against every unique element.
A = bsxfun(@eq,colvec(x),rowvec(UUID));

% Find which unique elements (columns of A) have a sum greater than 1
% (i.e., appear more than once in the original vector x).
take = (sum(A) > 1);
V = UUID(take);
% I retains only the columns corresponding to the duplicate values (V)
I = A(:,take);
end %fun