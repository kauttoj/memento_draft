function data = run_segment_preprocessor(input_nifti,input_motion,savepath,CFG)

%% General settings
cfg.overwrite = 0; % set to one if you are re-running preprocessing and you want to overwrite existing files
%cfg.bramilapath = '/triton/becs/scratch/braindata/shared/toolboxes/bramila/bramila'; % bramila toolbox path
cfg.bramilapath = '/triton/becs/scratch/braindata/kauttoj2/code/bramila_git/latest_bramila'; % bramila toolbox path
cfg.StdTemplate='/triton/becs/scratch/braindata/shared/HarvardOxford/MNI152_T1_2mm_brain.nii'; % 2mm MNI template from FSL
cfg.TR = CFG.TR; % TR from scanning protocol, used in bramila

%% Smoothing
cfg.do_spatial_smooth = 1;
cfg.smooth_FWHM = CFG.FWHM; % used in susan smoothing
cfg.smooth_method = 'FSLgauss'; % 'SPM', 'AFNI', 'FSL' or 'none'

cfg.FSLDIR = '';

cfg.max_tissue_pca_count = 0;

cfg.write = 0;    % write all intermediate EPI's
cfg.do_temporal_filtering = 1;

cfg.mot_derivatives = 0;  % motion regressor derivatives
cfg.detrend_type=CFG.DETREND_type;   % detrending type
cfg.filtertype = 'SPMhp';   % temporal filter type
cfg.filter_limits = [0,CFG.HP_limit,1,1];   % filter limits in Hz

cfg.infile =      input_nifti;
cfg.motionparam = input_motion;
cfg.motion_reg_type = 'standard';
cfg.save_path =   savepath;

%% RUN PREPROCESSING
cfg = bramila_checkparameters(cfg);

fprintf('\nFileID ''%s''\n',cfg.fileID);
fprintf('EPI path: %s\n',cfg.infile);
fprintf('Save folder: %s\n',cfg.outpath);

% create EPI mask based on signal quality
[cfg.vol,cfg.mask] = bramila_makemask(cfg);

% nullify all voxels outside the EPI mask
cfg.vol = bramila_maskdata(cfg);

% remove trends from the data
cfg.vol = bramila_detrend(cfg);

% regress out motion inside EPI mask
cfg.infile=[];
cfg = bramila_motionregress(cfg);

% temporal filtering
cfg = bramila_filter(cfg);

% spatial smoothing
cfg = bramila_smooth(cfg);

data = cfg.vol;

end