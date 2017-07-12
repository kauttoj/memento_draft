function compute_roi_correlations_stats_stage3(ITERATIONS,nworker)

if nargin<2
    nworker = 4;
end

diary('ROI_stats_diary.txt')

%S = {'S7','S8', 'S9'};%,'S10','S12','S13','S15','S16','S17','S19','S21','S22','S23'}; % S13 incomplete

MAX_IMAGES = [1369,1369,1369];

HOME = pwd;

TR = 1.56;
first_volume_offset = 0.15 + TR/2;

addpath('/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc');
TIMING_FILE = '/triton/becs/scratch/braindata/kauttoj2/Memento/2015/analysis/memento_git_project/Misc/chrono_original_segments/memento_chrono_original_color_timings_FINAL.txt';
COLOR_SEGMENTS = parse_timing_file(TIMING_FILE);

SELECTED_IDs=[2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,20,21,22];

ind=[];
for i=1:length(COLOR_SEGMENTS)
    if ismember(COLOR_SEGMENTS(i).ID,SELECTED_IDs)>0
        ind(end+1)=i;
    end
end
COLOR_SEGMENTS = COLOR_SEGMENTS(ind);
clear ind;

fprintf('\n-- Processing timing data ---\n')

averaged_volumes = struct();
k=0;
for i=1:length(COLOR_SEGMENTS)
    
    k=k+1;
    [fmri_volumes,volume_times,BOLD_envelope_fmri] = give_BOLD_envelope(COLOR_SEGMENTS(i).ORIGINAL_time(1),diff(COLOR_SEGMENTS(i).ORIGINAL_time),TR,first_volume_offset,0.5);
    
    fmri_volumes(1:2)=[];
    fmri_volumes(end-7:end)=[];
    
    averaged_volumes(k).volumes = fmri_volumes;
    averaged_volumes(k).session = COLOR_SEGMENTS(i).ORIGINAL_ses;
    averaged_volumes(k).id=COLOR_SEGMENTS(i).ID;
    averaged_volumes(k).type = 'color_segment';
    ind = averaged_volumes(k).volumes<=MAX_IMAGES(averaged_volumes(k).session);
    averaged_volumes(k).volumes=averaged_volumes(k).volumes(ind);
    
    averaged_volumes_size(k)=length(averaged_volumes(k).volumes);
    
    sessiondata(i)=averaged_volumes(k).session;
    
end

cd(HOME);

load('volume_selection_data.mat');

vols = 0;
for i=1:length(averaged_volumes_all)
    v=0;
    for j=1:length(averaged_volumes_all{i})
        v = v + length(averaged_volumes_all{i}(j).volumes);
    end
    v = v/length(averaged_volumes_all{i});
    vols = vols + v;
end
N_nullvols = round(vols/length(averaged_volumes_all));

ref_ind = round(length(averaged_volumes_all)/2);
N_keyframes = length(averaged_volumes_all{ref_ind});
for i=1:N_keyframes
    null_selection_session(i) = averaged_volumes_all{ref_ind}(i).session;
end

for i=1:3
   sesinds{i}=find(sessiondata==i);
   nses(i)=nnz(null_selection_session==i);
end

for iter = 1:ITERATIONS
    nulldata_selection=[];
    for ses=1:3
        r = randsample(length(sesinds{ses}),nnz(null_selection_session==ses),false)';
        ind = sesinds{ses}(r);      
        nulldata_selection=[nulldata_selection,ind];
    end
    
    nulldata(iter).startind=zeros(1,N_keyframes);
    nulldata(iter).session=zeros(1,N_keyframes);
    nulldata(iter).segment=zeros(1,N_keyframes);
    for i=1:N_keyframes
        r = randi(averaged_volumes_size(nulldata_selection(i))-N_nullvols+1);
        nulldata(iter).startind(i) = averaged_volumes(nulldata_selection(i)).volumes(r);
        nulldata(iter).session(i) = averaged_volumes(nulldata_selection(i)).session;
        nulldata(iter).segment(i) = averaged_volumes(nulldata_selection(i)).id;
    end
end

save('volume_selection_data_stats.mat','COLOR_SEGMENTS','nulldata','ref_ind','N_nullvols','averaged_volumes','ITERATIONS','-v7.3');

load('roi_and_mask_data.mat');

N_rois = length(my_masks);
for i=1:N_rois
    mask_voxels(i)=length(my_masks{i});
end

fprintf('\n--- Number of ROIs is %i ---\n',N_rois);

triu_inds = triu(ones(N_keyframes,N_keyframes),1)>0;

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

fprintf('\n--- Starting iterations ---\n');

parfor subj_ind = 1:length(S);
    
    subj_nullvals=nan(N_rois,ITERATIONS);
    
    subj = S{subj_ind};
    
    fprintf('\nLoading data for Subject %s... ',subj);
    
    cd(HOME);
    data = load_data([subj,'_roi_data.mat']);
    
    if length(data{1})~=N_rois
        error('ROI count does not match');
    end
    
    avg_img=cell(1,N_rois);
    for i=1:N_rois
        avg_img{i}=nan(mask_voxels(i),N_keyframes);
    end
    
    fprintf('done\nStarting iterations\n');    
            
    for iter = 1:ITERATIONS
        
        if mod(iter,500)==0
            fprintf('... iter %i/%i\n',iter,ITERATIONS);
        end
        BAD_FRAMES = [];
        
        for i = 1:N_keyframes
            ses = nulldata(iter).session(i);
            ind = nulldata(iter).startind(i) + (1:N_nullvols) - 1;
            
            MAX_TIMEPOINTS = size(data{ses}{1},1);
            
            if ind(end)>MAX_TIMEPOINTS
                BAD_FRAMES(end+1)=i;
            else
                for j=1:N_rois
                    avg_img{j}(:,i)=sum(data{ses}{j}(ind,:));
                end
            end
        end
        
        if ~isempty(BAD_FRAMES)
            iter
            warning('Bad frames found!');
        end
        
        for j=1:N_rois
            avg_img{j}(:,BAD_FRAMES)=[];
            avg_img{j}=avg_img{j}/N_nullvols;
            cormat = corr(avg_img{j});
            subj_nullvals(j,iter)=mean(cormat(triu_inds));
        end                                
        
    end
    
    if nnz(isnan(subj_nullvals))>0
        error('nullvals found!')
    end
    
    fprintf('done, saving data.\n');
    save_data([subj,'_nullvals.mat'],subj_nullvals);

end

subj_nullvals=[];

all_nullvals = zeros(N_rois,ITERATIONS);
for subj_ind = 1:length(S);
    
    subj = S{subj_ind};
    load([subj,'_nullvals.mat']);
    all_nullvals = all_nullvals+subj_nullvals;
    
end
all_nullvals=all_nullvals/length(S);

cd(HOME);
save('combined_nullvals','nulldata','all_nullvals','S','my_rois','my_masks','searchlightRad_mm','ITERATIONS','N_workers','-v7.3');

diary off

fprintf('\n--- ALL DONE!! ----\n')

function save_data(filename,subj_nullvals)

save(filename,'subj_nullvals','-v7.3');

function data = load_data(filename)

data=[];
load(filename);