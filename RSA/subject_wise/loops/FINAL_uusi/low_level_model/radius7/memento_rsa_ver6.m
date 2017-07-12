%% DEMO4_RSAsearchlight_sim
% simulates fMRI data for a number of subjects and runs searchlight
% analysis using RSA, computes the similarity maps for each subject
% and does group level inference.

% demo of the searchlight analysis
%%%%%%%%%%%%%%%%%%%%
%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%
clear all;
clc
close all;

returnHere = pwd; % We'll come back here later
toolboxRoot = '/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox';
addpath(genpath(toolboxRoot));

load('../../averaged_patterns_summary.mat');

for s=1:length(subj_averaged_volumes)
    
    nConditions = length(subj_averaged_volumes{s});
    label = [];
    a=ones(nConditions,nConditions);
    for i=1:nConditions
        label{i}=subj_averaged_volumes{s}(i).type;
        for j=1:nConditions
            if i==j
                a(i,j)=0;
            elseif ( (strcmp(subj_averaged_volumes{s}(i).type,'initial') && strcmp(subj_averaged_volumes{s}(j).type,'keyframe')) ...
                    || ( strcmp(subj_averaged_volumes{s}(i).type,'keyframe') && strcmp(subj_averaged_volumes{s}(j).type,'initial') ) ) ...
                    && subj_averaged_volumes{s}(i).id == subj_averaged_volumes{s}(j).id                
                a(i,j)=0;
            end
        end
    end
    all_models{s}(1).name = 'low-level_model';
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
userOptions.random_iterations = 1e+6;
userOptions.random_seed = 1;
userOptions.rootPath = pwd;
userOptions.betaPath = '/triton/becs/scratch/braindata/kauttoj2/code/RSAtoolbox/Demos/demoTest/[[subjectName]]/[[betaIdentifier]]';
userOptions.subjectNames = S;
%[userOptions.alternativeConditionLabels{1:nConditions}] = deal(' ');
%userOptions.useAlternativeConditionLabels = true;
%userOptions.conditionColours = rand(nConditions,3);
%userOptions.convexHulls = 1:nConditions;
%userOptions.distance = 'Correlation'; % NO EFFECT IN searchlightMapping_fMRI
%userOptions.distanceMeasure = 'Spearman'; % NO EFFECT IN searchlightMapping_fMRI
%userOptions.significanceTestPermutations = 1000;
userOptions.rankTransform = 1;
%                     userOptions.colourScheme = [64x3 double]
% userOptions.rubberbands = 1;
% userOptions.criterion = 'metricstress';
% userOptions.displayFigures = 1;
% userOptions.saveFiguresPDF = 1;
% userOptions.saveFiguresFig = 0;
% userOptions.saveFiguresPS = 0;
% userOptions.saveFiguresEps = 0;
%userOptions.nResamplings = 1000;
%userOptions.resampleSubjects = 1;
%userOptions.resampleConditions = 1;
userOptions.RDMname = 'referenceRDM';
%userOptions.plottingStyle = 2;
userOptions.voxelSize = [2,2,2];
userOptions.searchlightRadius = 7;

searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
%searchlightOptions.nConditions = nConditions; % given by the number of pattern volumes

searchlightOptions.searchlight_sphere = sphericalRelativeRoi(userOptions.searchlightRadius,userOptions.voxelSize);
searchlightOptions.N_searchlight = size(searchlightOptions.searchlight_sphere,1);

searchlightOptions.N_MIN_VOXELS = max(10,ceil(searchlightOptions.N_searchlight*0.5));

%load([returnHere,filesep,'sampleMask_org.mat'])
%load([returnHere,filesep,'anatomy.mat']);% load the resliced structural image
%models = constructModelRDMs(modelRDMs_SL_sim, userOptions);
%nCond = searchlightOptions.nConditions;
Nsubjects = length(userOptions.subjectNames);
if Nsubjects~=length(subj_averaged_volumes)
    error('number of subject does not match')
end

mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];
%% simulate the data and compute the correlation maps per subject

maskName = 'mask';

nii = load_nii('../../PATTERN_mask.nii');
mask = nii.img;

rs_all=cell(Nsubjects,1);
rs_all_res=cell(Nsubjects,1);
ps_all=cell(Nsubjects,1);
ns_all=cell(Nsubjects,1);
kloss_all=cell(Nsubjects,1);
null_rs_all=cell(Nsubjects,1);

parfor subI = 1:Nsubjects
    subject = ['subject',num2str(subI)];
    
    nii = load_nii(['../../',userOptions.subjectNames{subI},'_averaged_patterns.nii']);
    singleSubjectVols = nii.img;
    siz=size(singleSubjectVols);
    
    model = all_models{subI};
    
    if siz(4) ~= size(model.RDM,1)
        error('Wrong volume number!')
    end
    
    singleSubjectVols=reshape(singleSubjectVols,[],siz(4));
    % rs        4D array of 3D maps (x by y by z by model index) of
    %               correlations between the searchlight pattern similarity
    %               matrix and each of the model similarity matrices.
    %
    % ps        4D array of 3D maps (x by y by z by model index) of p
    %               values computed for each corresponding entry of smm_rs.
    %
    % ns             an array of the same dimensions as the volume, which
    %               indicates for each position how many voxels contributed
    %               data to the corresponding values of the infomaps.
    %               this is the number of searchlight voxels, except at the
    %               fringes, where the searchlight may illuminate voxels
    %               outside the input-data mask or voxel with all-zero
    %               time-courses (as can arise from head-motion correction).
    fprintf(['computing correlation maps for subject %i \n'],subI)    
    
    [rs_all{subI}, ps_all{subI}, ns_all{subI},null_rs_all{subI},rs_all_res{subI}] = searchlightMapping_fMRI(singleSubjectVols, model, mask, userOptions, searchlightOptions);
    %[rs_all{subI}, ps_all{subI}, ns_all{subI},kloss_all{subI},null_rs_all{subI}] = searchlightMapping_fMRI_extended(singleSubjectVols, models, mask, userOptions, searchlightOptions);
   
    nullcount(subI)=size(null_rs_all{subI},1);
    
    %rs_all{subI} = simple_glm_test(singleSubjectVols, models, mask, userOptions);
    
end

%save('kloss_all.mat','kloss_all')

gotoDir(userOptions.rootPath, 'Maps');
for subI = 1:Nsubjects
    subject = ['subject',num2str(subI)];
    
    rs=rs_all{subI};
    save(['rs_',subject,'.mat'],'rs');
    save_nii_oma(rs,['rs_',subject,'.nii']);
    
    rs_res=rs_all_res{subI};
    save(['rs_res_',subject,'.mat'],'rs_res');
    save_nii_oma(rs_res,['rs_res_',subject,'.nii']);
    
    ns=ns_all{subI};
    save(['ns_',subject,'.mat'],'ns');
    
    null_rs = null_rs_all{subI};
    save(['null_rs_',subject,'.mat'],'null_rs');
end

cd(returnHere)

%subplot(324);
% figure;
% image(scale01(rankTransform_equalsStayEqual(models(1).RDM,1)),'CDataMapping','scaled','AlphaData',~isnan(models(1).RDM));
% axis square off
% colormap(RDMcolormap)
% title('\bftested model RDM')


%% load the previously computed rMaps and concatenate across subjects
% prepare the rMaps:
supermask = mask;

for subjectI = 1:Nsubjects
    fprintf(['loading the correlation maps for subject %d \n'],subjectI);
    
    load([userOptions.rootPath,filesep,'Maps',filesep,'rs_subject',num2str(subjectI),'.mat'])
    rMaps{subjectI} = rs;
    
    load([userOptions.rootPath,filesep,'Maps',filesep,'ns_subject',num2str(subjectI),'.mat'])
    supermask=supermask.*(ns>=searchlightOptions.N_MIN_VOXELS);
        
    load([userOptions.rootPath,filesep,'Maps',filesep,'null_rs_subject',num2str(subjectI),'.mat'])  
    nullvalues{subjectI} = null_rs(1:min(nullcount),:);   
    
end

fprintf('\nSupermask size is %i (from %i)\n\n',nnz(supermask),nnz(mask));

supermask_ind = find(supermask);

% concatenate across subjects
modelI=1;

nullvalues_matrix = zeros(min(nullcount),Nsubjects);
%for modelI = 1:numel(models)
for subI = 1:Nsubjects
    thisRs = rMaps{subI};
    if subI==1
        thisModelSims=zeros([size(thisRs(:,:,:,modelI)),Nsubjects]);
    end
    thisModelSims(:,:,:,subI) = thisRs(:,:,:,modelI);
    
    save_nii_oma(squeeze(thisRs(:,:,:,modelI)),sprintf('S%i_corrmap.nii',subI));
    
    if min(nullcount)>1000
        nullvalues_matrix(:,subI)=nullvalues{subI}(:,1);
        
        a=thisRs(supermask_ind)';
        [Th_kauppi,Th_info_kauppi]=compute_pvals_kauppi(nullvalues_matrix(:,subI)', a);
        mean_map_kauppi = thisRs;
        mean_map_kauppi(mean_map_kauppi<Th_kauppi(2))=nan;
        mean_map_kauppi(isnan(mean_map_kauppi))=0;
        fprintf('Subject %i had %i voxels over threshold\n',subI,nnz(mean_map_kauppi));
        mean_map_kauppi= extentThreshold(mean_map_kauppi,50);
        save_nii_oma(mean_map_kauppi,sprintf('S%i_kauppi_p005_extent50.nii',subI));
    end
end
% obtain a pMaps from applying a 1-sided signrank test and also t-test to
% the model similarities:
nullvalues_mean = mean(nullvalues_matrix,2);

siz = size(supermask);
p1=nan(siz);
p2=nan(siz);
mean_map = nan(siz);

fprintf('\nComputing second-level statistics...')

%     ind = find(supermask>0)';
%     tempdata = reshape(thisModelSims,[],size(thisModelSims,4));
%     for i=ind
%         [~,p1(i)] = ttest(tempdata(i,:),0,0.05,'right');
%         p2(i) = signrank_onesided(tempdata(i,:));
%     end
%     mean_map(ind)=mean(tempdata(ind,:));
%
%     error('breakpoint here')
%

for x=1:siz(1)
    if mod(x,10)==0
        fprintf('%i ',x);
    end
    for y=1:siz(2)
        for z=1:siz(3)
            if supermask(x,y,z) > 0
                a = squeeze(thisModelSims(x,y,z,:));
                [~,p1(x,y,z)] = ttest(a,0,0.05,'right');
                p2(x,y,z) = signrank_onesided(a);
                mean_map(x,y,z)=mean(a);
            end
        end
    end
end
fprintf(' done\n')

% apply FDR correction
pThrsh_t = FDRthreshold(p1,0.05,supermask);
pThrsh_sr = FDRthreshold(p2,0.05,supermask);
p_bnf = 0.05/nnz(supermask);
% mark the suprathreshold voxels in yellow

supraThreshMarked_t = zeros(size(p1));
ind_t = (p1 <= pThrsh_t);
supraThreshMarked_t(ind_t) = 1;
meanmap_t = zeros(size(p1));
meanmap_t(ind_t)=mean_map(ind_t);
meanmap_t_p0001 = zeros(size(p1));
meanmap_t_p0001(p1 < 0.001)=mean_map(p1 < 0.001);

supraThreshMarked_sr = zeros(size(p2));
ind_sr = (p2 <= pThrsh_sr);
supraThreshMarked_sr(ind_sr) = 1;
meanmap_sr = zeros(size(p2));
meanmap_sr(ind_sr)=mean_map(ind_sr);
meanmap_sr_p0001 = zeros(size(p2));
meanmap_sr_p0001(p2 < 0.001)=mean_map(p2 < 0.001);

if min(nullcount)>1000
    [Th_kauppi,Th_info_kauppi]=compute_pvals_kauppi(nullvalues_mean',mean_map(supermask_ind)');
    mean_map_kauppi = mean_map;
    mean_map_kauppi(mean_map_kauppi<Th_kauppi(2))=nan;
    mean_map_kauppi(isnan(mean_map_kauppi))=0;
    mean_map_kauppi= extentThreshold(mean_map_kauppi,50);
    save_nii_oma(mean_map_kauppi,'meanmap_kauppi_p005_extent50.nii');
end
%     % display the location where the effect was inserted (in green):
%     brainVol = addRoiToVol(map2vol(anatVol),mask2roi(mask),[1 0 0],2);
%     brainVol_effectLoc = addBinaryMapToVol(brainVol,mask,[0 1 0]);
%     showVol(brainVol_effectLoc,'simulated effect [green]',2);
%     handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO4_simulatedEffectRegion'],userOptions);
%
%     % display the FDR-thresholded maps on a sample anatomy (signed rank test) :
%     brainVol = addRoiToVol(map2vol(anatVol),mask2roi(mask),[1 0 0],2);
%     brainVol_sr = addBinaryMapToVol(brainVol,supraThreshMarked_sr.*mask,[1 1 0]);
%     showVol(brainVol_sr,'signrank, E(FDR) < .05',3)
%     handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO4_signRank'],userOptions);
%
%     % display the FDR-thresholded maps on a sample anatomy (t-test) :
%     brainVol = addRoiToVol(map2vol(anatVol),mask2roi(mask),[1 0 0],2);
%     brainVol_t = addBinaryMapToVol(brainVol,supraThreshMarked_t.*mask,[1 1 0]);
%     showVol(brainVol_t,'t-test, E(FDR) < .05',4)
%     handleCurrentFigure([returnHere,filesep,'DEMO4',filesep,'results_DEMO2_tTest'],userOptions);


supraThreshMarked_t=extentThreshold(supraThreshMarked_t,50);
supraThreshMarked_sr=extentThreshold(supraThreshMarked_sr,50);

meanmap_t = extentThreshold(meanmap_t,50);
meanmap_sr = extentThreshold(meanmap_sr,50);

meanmap_t_p0001 = extentThreshold(meanmap_t_p0001,100);
meanmap_sr_p0001 = extentThreshold(meanmap_sr_p0001,100);

save_nii_oma(supraThreshMarked_t,'supraThreshMarked_t_p005_extent50.nii');
save_nii_oma(supraThreshMarked_sr,'supraThreshMarked_sr_p005_extent50.nii');

save_nii_oma(meanmap_t,'meanmap_t_p005_extent50.nii');
save_nii_oma(meanmap_sr,'meanmap_sr_p005_extent50.nii');

save_nii_oma(meanmap_t_p0001,'meanmap_t_p0001noncor_extent100.nii');
save_nii_oma(meanmap_sr_p0001,'meanmap_sr_p0001noncor_extent100.nii');

save_nii_oma(mean_map,'meanmap.nii');

save_nii_oma(supermask,'supermask.nii');

%end
cd(returnHere);

save('memento_rsa_results.mat','pThrsh_t','pThrsh_sr','p_bnf','supraThreshMarked_t','supraThreshMarked_sr','p1','p2','meanmap_t','meanmap_sr','mean_map','searchlightOptions','userOptions','mask','supermask','nullvalues_mean','Th_kauppi','Th_info_kauppi','all_models','S');

disp('---- ALL DONE !!! ------')
