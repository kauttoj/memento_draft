function mask = compute_analysis_mask(data_path,S)
%COMPUTE_ANALYSIS_MASK Summary of this function goes here
%   Detailed explanation goes here

HOME = pwd;

mask = ones(91,109,91);

for ses=1:length(data_path)
    cd(data_path{ses});
    for s=S
        
        nii=load_nii([s{1},'_analysis_mask.nii']);
        mask = mask.*nii.img;
        
    end
end

cd(HOME);

end

