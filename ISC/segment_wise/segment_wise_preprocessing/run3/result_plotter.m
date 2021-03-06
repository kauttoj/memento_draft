clc;

N = length(ISC_matrices);

for i=1:N
        
    str = sprintf('segment_ID_%i',MOVIE_SEGMENTS(i).ID);    
    
    ind = ISC_perspective_pval_cor{i}<0.05;
    val = ISC_perspective_corr{i}(ind);    
    map = 0*group_mask;    
    map(group_mask_ind(ind))=val;
    map=extentThreshold(map,20);
    if nnz(map)>0
        save_nii_oma(map,[str,'_perspective_isc_diff.nii']);
        %save_nii_oma(-map.*(map<0),[str,'_perspective_isc_diff_NEG.nii']);
    end
    
    N1 = nnz(map);
    
    ind = ISC_perspective_pval{i}<0.005;
    val = ISC_perspective_corr{i}(ind);
    map = 0*group_mask;    
    map(group_mask_ind(ind))=val;
    map=extentThreshold(map,50);
    if nnz(map)>0
        save_nii_oma(map,[str,'_perspective_isc_diff_uncorrected.nii']);
        %save_nii_oma(-map.*(map<0),[str,'_perspective_isc_diff_NEG.nii']);
    end
        
    ind = ISC_difference_stats{i}.twotail_pvals_cor<0.05;
    val = ISC_difference_stats{i}.tvals(ind);
    map = 0*group_mask;    
    map(group_mask_ind(ind))=val;
    map=extentThreshold(map,20);
    
    N2 = nnz(map);
    
    tempmap = map.*(map>0);
    if nnz(tempmap)>0
        save_nii_oma(tempmap,[str,'_isc_diff_POS.nii']);
    end
    tempmap = map.*(map<0);
    if nnz(tempmap)>0
        save_nii_oma(-tempmap,[str,'_isc_diff_NEG.nii']);    
    end    
    
    map = 0*group_mask;   
    map(group_mask_ind) = KRON_ISC_vals;
    save_nii_oma(map,[str,'_isc_vals_KRON.nii']);
    
    map = 0*group_mask;   
    map(group_mask_ind) = ORIGINAL_ISC_vals;
    save_nii_oma(map,[str,'_isc_vals_ORIGINAL.nii']);
            
    fprintf('\n>>>>>>>> ID %i: %i voxels for DIFFERENCE, %i voxels for PERSPECTIVE\n\n',MOVIE_SEGMENTS(i).ID,N2,N1);
    
end