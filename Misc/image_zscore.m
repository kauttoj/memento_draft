function mat = image_zscore(mat)
%IMAGE_ZSCORE Summary of this function goes here
%   Detailed explanation goes here

mat = mat - mean2(mat);
mat = mat/std(mat(:));

end

