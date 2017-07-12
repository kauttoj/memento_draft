
    nii_p=load_nii('NEG_OneSampT_tfce_corrp_tstat1.nii');
    nii_t=load_nii('NEG_OneSampT_tstat1.nii');    
    img = nii_t.img.*(nii_p.img>0.975); 
    save_nii_oma(img,'GLM_image_NEG.nii');
    files{3}=[pwd,filesep,'GLM_image_NEG.nii'];
    
    nii_p=load_nii('POS_OneSampT_tfce_corrp_tstat1.nii');
    nii_t=load_nii('POS_OneSampT_tstat1.nii');    
    img = nii_t.img.*(nii_p.img>0.975); 
    save_nii_oma(img,'GLM_image_POS.nii');
    files{2}=[pwd,filesep,'GLM_image_POS.nii']; 
    
    nii_p=load_nii('OneSampT_tfce_corrp_tstat1.nii');
    nii_t=load_nii('OneSampT_tstat1.nii');    
    img = nii_t.img.*(nii_p.img>0.99); 
    save_nii_oma(img,'RSA_image_POS.nii');
    files{1}=[pwd,filesep,'RSA_image_POS.nii'];    
    
    filename = create_mricroGL_script('RSA_GLM_combined',files,pwd);