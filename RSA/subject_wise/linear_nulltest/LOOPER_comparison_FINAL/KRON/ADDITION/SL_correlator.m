function rs_all = SL_correlator(FOLDER,RADIUS)

cd(FOLDER);
returnHere = pwd; % We'll come back here later
toolboxRoot = '/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox';
addpath(genpath(toolboxRoot));

load('averaged_patterns_summary.mat');

for s=1:length(subj_averaged_volumes)
    
    nConditions = length(subj_averaged_volumes{s});
    label = [];
    a=ones(nConditions,nConditions);
    for i=1:nConditions
        label{i}=subj_averaged_volumes{s}(i).type;
        for j=1:nConditions
            if i==j
                a(i,j)=0;
            elseif strcmp(subj_averaged_volumes{s}(i).type,'keyframe') && strcmp(subj_averaged_volumes{s}(j).type,'keyframe')
                a(i,j)=0;
            end
        end
    end
    
    if nnz(a~=0)>0
        error('Incorrect model!')
    end
    
    all_models{s}(1).name = 'key-frame_model';
    all_models{s}(1).RDM = a;
    all_models{s}(1).label = label;    
    
end

nworker = 7;

% create local cluster
try
    myCluster = gcp('nocreate');
    if isempty(myCluster)
        %delete(gcp)
        myCluster = parcluster('local');
        myCluster.NumWorkers=nworker;
        parpool(myCluster);
    end
    N_workers = myCluster.NumWorkers;
catch err % old matlab?
    if ~matlabpool('size')
        eval(['matlabpool local ',num2str(nworker)]);
    end
    N_workers = matlabpool('size');
end

userOptions.RDM_correlation_type='Spearman';
userOptions.seed = 666;
userOptions.analysisName = 'mementoRSA';
userOptions.random_iterations = 100000;
userOptions.random_seed = 1;
userOptions.rootPath = pwd;
userOptions.betaPath = '/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox/Demos/demoTest/[[subjectName]]/[[betaIdentifier]]';
userOptions.subjectNames = S;
userOptions.rankTransform = 1;
userOptions.RDMname = 'referenceRDM';
userOptions.voxelSize = [2,2,2];
userOptions.searchlightRadius = RADIUS;

searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;

searchlightOptions.searchlight_sphere = sphericalRelativeRoi(userOptions.searchlightRadius,userOptions.voxelSize);
searchlightOptions.N_searchlight = size(searchlightOptions.searchlight_sphere,1);

Nsubjects = length(userOptions.subjectNames);
if Nsubjects~=length(subj_averaged_volumes)
    error('number of subject does not match')
end

mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];
%% simulate the data and compute the correlation maps per subject

maskName = 'mask';

nii = load_nii('PATTERN_mask.nii');
mask = nii.img;

rs_all=cell(Nsubjects,1);
rs_all_res=cell(Nsubjects,1);
ps_all=cell(Nsubjects,1);
ns_all=cell(Nsubjects,1);
null_rs_all=cell(Nsubjects,1);

parfor subI = 1:Nsubjects
    subject = ['subject',num2str(subI)];
    
    nii = load_nii([userOptions.subjectNames{subI},'_averaged_patterns.nii']);
    singleSubjectVols = nii.img;
    siz=size(singleSubjectVols);
    
    model = all_models{subI};
    
    if siz(4) ~= size(model.RDM,1)
        error('Wrong volume number!')
    end
    
    singleSubjectVols=reshape(singleSubjectVols,[],siz(4));

    fprintf(['computing correlation maps for subject %i \n'],subI)
        
    [rs_all{subI},ps_all{subI}, ns_all{subI}] = SL_compute(singleSubjectVols, model, mask, userOptions, searchlightOptions);
    %[rs_all{subI}, ps_all{subI}, ns_all{subI},kloss_all{subI},null_rs_all{subI}] = searchlightMapping_fMRI_extended(singleSubjectVols, models, mask, userOptions, searchlightOptions);
       
    %rs_all{subI} = simple_glm_test(singleSubjectVols, models, mask, userOptions);
    
end


%save('kloss_all.mat','kloss_all')

gotoDir(userOptions.rootPath, 'Maps');
for subI = 1:Nsubjects
    subject = ['subject',num2str(subI)];
    
    rs=rs_all{subI};
    save(['rs_',subject,'.mat'],'rs');
    save_nii_oma(rs,['rs_',subject,'.nii']);
    
    ns=ns_all{subI};
    save(['ns_',subject,'.mat'],'ns');
end

cd(returnHere)

end

function [smm_rs, smm_ps, n, null_rs,smm_rs_res] = SL_compute(fullBrainVolumes, models, mask, userOptions, localOptions)
% THIS IS THE LATEST VERSION!
% ARGUMENTS
% fullBrainVolumes	A voxel x condition x session matrix of activity
% 				patterns.
%
% models		A struct of model RDMs.
%
% mask     		A 3d or 4d mask to perform the searchlight in.
%
% userOptions and localOptions
%
% RETURN VALUES
% smm_rs        4D array of 3D maps (x by y by z by model index) of
%               correlations between the searchlight pattern similarity
%               matrix and each of the model similarity matrices.
%
% smm_ps        4D array of 3D maps (x by y by z by model index) of p
%               values computed for each corresponding entry of smm_rs.
%
% n             an array of the same dimensions as the volume, which
%               indicates for each position how many voxels contributed
%               data to the corresponding values of the infomaps.
%               this is the number of searchlight voxels, except at the
%               fringes, where the searchlight may illuminate voxels
%               outside the input-data mask or voxel with all-zero
%               time-courses (as can arise from head-motion correction).
%
% mappingMask_actual
%               3D mask indicating locations for which valid searchlight
%               statistics have been computed.
%
% Based on Niko Kriegeskorte's searchlightMapping_RDMs.m
%
% Additions by Cai Wingfield 2-2010:
% 	- Now skips points in the searchlight where there's only one voxel inside.
% 	- Now takes a userOptions struct for the input parameters.

localOptions = setIfUnset(localOptions, 'averageSessions', true);

if isfield(localOptions,'N_MIN_VOXELS')
    N_MIN_VOXELS = localOptions.N_MIN_VOXELS;
else
    N_MIN_VOXELS = 5;
end

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',userOptions.seed));

%% Figure out whether to average over sessions or not
if localOptions.averageSessions
    for sessionNumber = 1:size(fullBrainVolumes,3)
        thisSessionId = ['s' num2str(sessionNumber)];
        t_patsPerSession.(thisSessionId) = fullBrainVolumes(:,:,sessionNumber)';
    end%for:sessionNumber
else
    justThisSession = 1;
    t_pats = fullBrainVolumes(:,:,justThisSession)';
    
    fprintf(['\nYou have selected not to average over sessions.\n         Only session number ' num2str(justThisSession) ' will be used.\n']);
    
end%if

if isfield(userOptions,'RDM_correlation_type')
    RDM_CORRELATION_TYPE = userOptions.RDM_correlation_type;
else
    RDM_CORRELATION_TYPE = 'Spearman';
end
fprintf('\nRDM correlation method: %s\n',RDM_CORRELATION_TYPE);

%% Get parameters
voxSize_mm = userOptions.voxelSize;
searchlightRad_mm = userOptions.searchlightRadius;
monitor = localOptions.monitor;
nConditions = size(fullBrainVolumes, 2);

clear fullBrainVolumes;

temp_models=models;
% Prepare models
for i=1:length(temp_models)
    siz=size(temp_models(i).RDM);
    a=ones(siz);
    ind = find(a);
    a(ind)=1:length(ind);
    temp_models(i).RDM=a;
    modelRDMs_ltv_MATRICES{i}=models(i).RDM;
end
modelRDMs_ltv_IND = permute(unwrapRDMs(vectorizeRDMs(temp_models)), [3 2 1]);
for i=1:size( modelRDMs_ltv_IND,2)
    if ~all(modelRDMs_ltv_IND(:,i)-modelRDMs_ltv_IND(1,i) == 0)
        error('matrix indices not equal!')
    end
end
modelRDMs_ltv_IND = modelRDMs_ltv_IND(1,:);

modelRDMs_ltv = permute(unwrapRDMs(vectorizeRDMs(models)), [3 2 1]);

% Prepare masks
mask(isnan(mask)) = 0; % Just in case!
if ndims(mask)==3
    inputDataMask=logical(mask);
    mappingMask_request=logical(mask);
else
    inputDataMask=logical(mask(:,:,:,1));
    mappingMask_request=logical(mask(:,:,:,2));
end

% Check to see if there's more data than mask...
if localOptions.averageSessions
    for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
        thisSessionId = ['s' num2str(sessionNumber)];
        t_patsPerSession.(thisSessionId) = t_patsPerSession.(thisSessionId)(:, inputDataMask(:));
    end%for:sessionNumber
else
    if (size(t_pats,2)>sum(inputDataMask(:)))
        t_pats=t_pats(:,inputDataMask(:));
    end%if
end%if

% Other data
volSize_vox=size(inputDataMask);
nModelRDMs=size(modelRDMs_ltv,1);
rad_vox=searchlightRad_mm./voxSize_mm;
minMargin_vox=floor(rad_vox);


%% create spherical multivariate searchlight
[x,y,z]=meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
sphere=((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(searchlightRad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
sphereSize_vox=[size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )

if monitor, figure(50); clf; showVoxObj(sphere); end % show searchlight in 3D

% compute center-relative sphere SUBindices
[sphereSUBx,sphereSUBy,sphereSUBz]=ind2sub(sphereSize_vox,find(sphere)); % (SUB)indices pointing to sphere voxels
sphereSUBs=[sphereSUBx,sphereSUBy,sphereSUBz];
ctrSUB=sphereSize_vox/2+[.5 .5 .5]; % (c)en(t)e(r) position (sphere necessarily has odd number of voxels in each dimension)
ctrRelSphereSUBs=sphereSUBs-ones(size(sphereSUBs,1),1)*ctrSUB; % (c)en(t)e(r)-relative sphere-voxel (SUB)indices

nSearchlightVox=size(sphereSUBs,1);


%% define masks
validInputDataMask=inputDataMask;

if localOptions.averageSessions
    for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
        thisSessionId = ['s' num2str(sessionNumber)];
        sumAbsY=sum(abs(t_patsPerSession.(thisSessionId)),1);
    end%for:sessionNumber
else
    sumAbsY=sum(abs(t_pats),1);
end%if

validYspace_logical= (sumAbsY~=0) & ~isnan(sumAbsY); clear sumAbsY;
validInputDataMask(inputDataMask)=validYspace_logical; % define valid-input-data brain mask

if localOptions.averageSessions
    for sessionNumber = 1:numel(fieldnames(t_patsPerSession))
        thisSessionId = ['s' num2str(sessionNumber)];
        t_patsPerSession.(thisSessionId) = t_patsPerSession.(thisSessionId)(:,validYspace_logical);
        nVox_validInputData=size(t_patsPerSession.(thisSessionId),2);
    end%for:sessionNumber
else
    t_pats=t_pats(:,validYspace_logical); % reduce t_pats to the valid-input-data brain mask
    nVox_validInputData=size(t_pats,2);
end%if

mappingMask_request_INDs=find(mappingMask_request);
nVox_mappingMask_request=length(mappingMask_request_INDs);

if monitor
    disp([num2str(round(nVox_mappingMask_request/prod(volSize_vox)*10000)/100),'% of the cuboid volume requested to be mapped.']);
    disp([num2str(round(nVox_validInputData/prod(volSize_vox)*10000)/100),'% of the cuboid volume to be used as input data.']);
    disp([num2str(nVox_validInputData),' of ',num2str(sum(inputDataMask(:))),' declared input-data voxels included in the analysis.']);
end

volIND2YspaceIND=nan(volSize_vox);
volIND2YspaceIND(validInputDataMask)=1:nVox_validInputData;

% n voxels contributing to infobased t at each location
n=nan(volSize_vox);

%% similarity-graph-map the volume with the searchlight
smm_bestModel=nan(volSize_vox);
smm_ps=nan([volSize_vox,nModelRDMs]);
smm_rs=nan([volSize_vox,nModelRDMs]);
smm_rs_res=nan([volSize_vox,nModelRDMs]);
%searchlightRDMs = nan([nConditions, nConditions, volSize_vox]);

if monitor
    h_progressMonitor=progressMonitor(1, nVox_mappingMask_request,  'Similarity-graph-mapping...');
end

%%

[xx,yy,zz]=ind2sub(volSize_vox,mappingMask_request_INDs);

n_voxelcount=nan(1,nVox_mappingMask_request);

for cMappingVoxI=1:nVox_mappingMask_request
       
    x=xx(cMappingVoxI);
    y=yy(cMappingVoxI);
    z=zz(cMappingVoxI);    
    
    % compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
    cIllVoxSUBs=repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1])+ctrRelSphereSUBs;
    
    % exclude out-of-volume voxels
    outOfVolIs=(cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
        cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
        cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));
    
    cIllVoxSUBs=cIllVoxSUBs(~outOfVolIs,:);
    
    % list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
    cIllVox_volINDs=sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));
    
    % restrict searchlight to voxels inside validDataBrainMask
    cIllValidVox_volINDs=cIllVox_volINDs(validInputDataMask(cIllVox_volINDs));
    cIllValidVox_YspaceINDs=volIND2YspaceIND(cIllValidVox_volINDs);
    
    % note how many voxels contributed to this locally multivariate stat
    n_voxelcount(cMappingVoxI)=length(cIllValidVox_YspaceINDs);
end

if nnz(isnan(n_voxelcount))>0
    error('NaN''s in voxelcount!')
end
%%

N_iter = userOptions.random_iterations;
null_rs = nan(N_iter,size(modelRDMs_ltv,1));

rng(userOptions.random_seed);

%% THE BIG LOOP! %%

[xx,yy,zz]=ind2sub(volSize_vox,mappingMask_request_INDs);

if nVox_mappingMask_request~=length(xx)
    error('Something went wrong!')
end   

fprintf('\n Starting searchlight analysis ')

for cMappingVoxI=1:nVox_mappingMask_request
    
    if mod(cMappingVoxI,5000)==0
        if monitor
            progressMonitor(cMappingVoxI, nVox_mappingMask_request, 'Searchlight mapping Mahalanobis distance...', h_progressMonitor);
            %                 cMappingVoxI/nVox_mappingMask_request
        else
            fprintf('.');
        end%if
    end%if
    
    x=xx(cMappingVoxI);
    y=yy(cMappingVoxI);
    z=zz(cMappingVoxI);    
    
    % compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
    cIllVoxSUBs=repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1])+ctrRelSphereSUBs;
    
    % exclude out-of-volume voxels
    outOfVolIs=(cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
        cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
        cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));
    
    cIllVoxSUBs=cIllVoxSUBs(~outOfVolIs,:);
    
    % list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
    cIllVox_volINDs=sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));
    
    % restrict searchlight to voxels inside validDataBrainMask
    cIllValidVox_volINDs=cIllVox_volINDs(validInputDataMask(cIllVox_volINDs));
    cIllValidVox_YspaceINDs=volIND2YspaceIND(cIllValidVox_volINDs);
    
    % note how many voxels contributed to this locally multivariate stat
    n(x,y,z)=length(cIllValidVox_YspaceINDs);
    
    if n(x,y,z) < 2, continue; end%if % This stops the function crashing if it accidentally encounters an out-of-brain floating voxel (these can occur if, for example, skull stripping fails)
    
    if localOptions.averageSessions
        searchlightRDM = zeros(nConditions,nConditions);
        for session = 1:localOptions.nSessions
            sessionId = ['s' num2str(session)];
            searchlightRDM = searchlightRDM + squareform(pdist(t_patsPerSession.(sessionId)(:,cIllValidVox_YspaceINDs),'correlation'));
        end%for:sessions
        searchlightRDM = searchlightRDM / localOptions.nSessions;
    else
        searchlightRDM = squareform(pdist(t_pats(:,cIllValidVox_YspaceINDs), 'correlation'));
    end%if
    
    searchlightRDM = 1 - vectorizeRDM(searchlightRDM);
    
    rs = mean(searchlightRDM);
    
    ps = nan;
    
    if localOptions.fisher
        for i = 1:numel(rs)
            rs(i) = fisherTransform(rs(i));
        end%for:i
    end%if
            
    smm_ps(x,y,z,:) = ps;
    smm_rs(x,y,z,:) = rs;

    
end%for:cMappingVoxI

%% END OF THE BIG LOOP! %%
fprintf('\n');


end%function
