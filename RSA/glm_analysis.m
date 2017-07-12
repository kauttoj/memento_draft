function [vals,pvals,pvals_corrected,resultmap,mask_out,corrmap_res]=glm_analysis(average_data,datapath,mask,TAIL,FWHM)
%GLM_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here

if nargin<5
    FWHM=6;
elseif nargin<4
    TAIL='both';
    FWHM=6;
end

if ~(strcmp(TAIL,'right') || strcmp(TAIL,'left') || strcmp(TAIL,'both'))
    warning('Incorrect tail parameter')
    TAIL = 'right';
end

bramilapath = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila';
addpath(bramilapath);

load(average_data);

if ischar(mask)
    nii=load_nii(mask);
    mask=nii.img;
end
mask_ind = find(mask);
mask_out = mask;

fprintf('GLM analysis: %i subjects, mask size is %i voxels\n',length(subj_averaged_volumes),length(mask_ind));

for s=1:length(subj_averaged_volumes)    
    nConditions = length(subj_averaged_volumes{s});
    label = [];
    a=zeros(1,nConditions);
    for i=1:nConditions
        label{i}=subj_averaged_volumes{s}(i).type;
        if strcmp(subj_averaged_volumes{s}(i).type,'keyframe')
           a(i)=1;           
        end
    end
    regressors{s}=a;
    labels{s}=label;
end

Nsubjects = length(S);
corrmap = nan(length(mask_ind),Nsubjects);

needs_smoothing = zeros(1,Nsubjects);

for subI = 1:Nsubjects
    subject = ['subject',num2str(subI)];        
    
    file = [datapath,filesep,S{subI},'_averaged_patterns_smoothed.nii'];
    
    fprintf('Processing subject %s\n',S{subI});
    
    if ~exist(file,'file')
        needs_smoothing(subI)=1;
        fprintf('...smoothing data\n');   
        clear cfg;
        cfg.bramilapath = bramilapath;
        cfg.FSLDIR='';
        cfg.infile = [datapath,filesep,S{subI},'_averaged_patterns.nii'];
        cfg.smooth_FWHM = FWHM;
        cfg.smooth_method='SPM';
        cfg.write=0;
        cfg.infile_orig=cfg.infile;
        cfg = bramila_smooth(cfg);
        singleSubjectVols = cfg.vol;
        save_nii_oma(singleSubjectVols,file);
    else    
        fprintf('...loading data\n'); 
        nii = load_nii(file);
        singleSubjectVols = nii.img;
    end
    
    fprintf('...computing correlations\n');     
    siz=size(singleSubjectVols);    
    singleSubjectVols=reshape(singleSubjectVols,[],siz(4))';
    singleSubjectVols = singleSubjectVols(:,mask_ind);
    corrmap(:,subI)=corr(singleSubjectVols,regressors{subI}');   
    
    if nnz(isnan(corrmap(:,subI)))>0
        subject
        error('NaN''s found!!')
    end
    
end

if length(unique(needs_smoothing))>1
    warning('Some of the data was already smoothed! All should have same method and FWHM!')
end

corrmap = atanh(corrmap)';
[~,pvals] = ttest(corrmap,0,'tail',TAIL);
vals = mean(corrmap);
pvals_corrected = mafdr(pvals,'BHFDR',true);
resultmap = 0*mask;
resultmap(mask_ind(pvals_corrected<0.05))=vals(pvals_corrected<0.05);

corrmap_res = nan([size(mask),Nsubjects]);
for subI = 1:Nsubjects
    a=0*mask;
    a(mask_ind)=corrmap(subI,:);
    corrmap_res(:,:,:,subI)=a;
end



