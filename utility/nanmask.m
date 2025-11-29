function img = nanmask(img,mask)
% NANMASK applies a logical mask to an image array (which may be 2D, 3D,
%   or 4D) and sets all elements outside of the mask (where mask is false)
%   to Not a Number (NaN). This effectively zeros out irrelevant regions
%   of the image data using the NaN flag, preserving the valid regions.
%
%   Inputs:
%       img (numeric array): The input image data (can be 2D, 3D, or 4D).
%       mask (logical array): A 2D logical mask (true for valid pixels,
%                             false for pixels to be set to NaN). The mask
%                             must have the same 2D dimensions (rows/cols) as img.
%
%   Output:
%       img (numeric array): The masked image array, with elements outside
%                            the mask set to NaN.
%
%   Copyright (c) 2025 Christian Paolo Richter

% Determine the dimensions of the input image 'img'.
% size(img) returns [dimY, dimX, dimZ, dimC, ...].
% We only care about the third (Z) and fourth (C) dimensions for replication.
% The first two (Y and X) are discarded (~) as they match the mask dimensions.
[~,~,dimZ,dimC] = size(img);

% Masking operation:
% 1. not(mask): Inverts the 2D logical mask. Pixels to be preserved are now false;
%    pixels to be set to NaN are now true.
% 2. repmat(..., 1, 1, dimZ, dimC): Replicates the inverted 2D mask along the
%    third (Z-slices/time) and fourth (C-channels/color) dimensions to match
%    the full dimensionality of 'img'. This creates a mask of the same size as 'img'.
% 3. img(repmat(...)) = nan;: Uses the replicated inverted mask for logical
%    indexing into 'img'. All elements in 'img' corresponding to 'true' in the
%    replicated mask (i.e., outside the original 'mask' region) are set to NaN.
img(repmat(not(mask),1,1,dimZ,dimC)) = nan;

end %fun