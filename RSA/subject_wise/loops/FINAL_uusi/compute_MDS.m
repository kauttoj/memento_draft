clc;
close all;
clear variables;

load averaged_patterns_summary;

nii=load_nii('PATTERN_mask.nii');
mask = nii.img;

coord = [...
    51 31 56;...
    72 36 56;...
    65 71 60;...
    ...%45 53 52;...
    49 72 61];

%coord = [65,39,65;...
%    68,66,58;...
%    44,51,51;...
%    52,34,54];
searchlightRad_mm = 6;
voxSize_mm=[2,2,2];

%%-----------------------
% Other data
rad_vox=searchlightRad_mm./voxSize_mm;
minMargin_vox=floor(rad_vox);
% create spherical multivariate searchlight
[x,y,z]=meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
sphere=((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(searchlightRad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
sphereSize_vox=[size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )
% compute center-relative sphere SUBindices
[sphereSUBx,sphereSUBy,sphereSUBz]=ind2sub(sphereSize_vox,find(sphere)); % (SUB)indices pointing to sphere voxels
sphereSUBs=[sphereSUBx,sphereSUBy,sphereSUBz];
ctrSUB=sphereSize_vox/2+[.5,.5,.5]; % (c)en(t)e(r) position (sphere necessarily has odd number of voxels in each dimension)
ctrRelSphereSUBs=sphereSUBs-ones(size(sphereSUBs,1),1)*ctrSUB; % (c)en(t)e(r)-relative sphere-voxel (SUB)indices
nSearchlightVox=size(sphereSUBs,1);
%%-----------------------

nROIs=size(coord,1);

for j=1:nROIs
    ROI = 0*mask;
    for i=1:nSearchlightVox
        ROI(ctrRelSphereSUBs(i,1)+coord(j,1),ctrRelSphereSUBs(i,2)+coord(j,2),ctrRelSphereSUBs(i,3)+coord(j,3))=1;
    end
    ROIs{j}=ROI.*mask;
    
    [result,roilabel{j}] = findClusterLabels_total(ROIs{j},1,1,1);
end

MDS_coord=cell(1,nROIs);
for i=1:length(MDS_coord)
   MDS_coord{i}=nan(30,2,1); 
   mean_D{i}=0;
end

k=0;
for i=1:length(S)
    [ts_data,PATTERN_mask_ind] = load_nii_mask([S{i},'_averaged_patterns.nii'],ROIs,500);
    
    if size(ts_data{1},1)==30
        k=k+1;
        pois = 1:length(subj_averaged_volumes{i});
        ind=[];
        for j=1:length(subj_averaged_volumes{i})
            if strcmp(subj_averaged_volumes{i}(j).type,'keyframe')
                ind(end+1)=j;
            end
        end
        pois(ind)=[];
        pois=[];
        for l=1:nROIs

            ts_data{l}(pois,:)=[];
            
            D = pdist(ts_data{l},'correlation');
            [Y,eigvals] = cmdscale(D);
            
            mean_D{l}=mean_D{l}+D;
            %scaled_eigs(1:2,i)=eigvals(1:2)/max(abs(eigvals));
            
            if size(Y,1)==30
                MDS_coord{l}(:,:,k) = [Y(:,1),Y(:,2)];
            else
                error('!!!');
            end
        end
        
    else
        
        disp('');
    end            
end

mean_D_orig = mean_D;

%cd procrustes
for l=1:nROIs
    mean_D{l}=mean_D{l}/k;  
    
    [mean_MDS_coord{l}, stress, disparities] = mdscale(mean_D{l},2,...
        'criterion','metricstress','options', struct('MaxIter', 1000));
    
    
%     [Y,eigvals] = cmdscale(mean_D{l});
%     mean_eigvals{l}=eigvals/max(eigvals);
%     mean_MDS_coord{l}=[Y(:,1),Y(:,2)];
%     [mean_shape{l},pc_shape,std_shape,pc_projection,new_pt_out]=tangent_pca_shape(MDS_coord{l});
end
%cd ..


coord_mni=ind2mni(coord);

figure;
for i=1:size(coord,1)
    subplot(2,2,i);
    hold on;
    title(sprintf('%s\n(x=%i,y=%i,z=%i)',roilabel{i}{1},coord_mni(i,1),coord_mni(i,2),coord_mni(i,3)));
    box on;
    for j=1:size(mean_MDS_coord{i},1)
        if mod(j,2)==0
            col='r';
        else
            col='b';
        end
        plot(mean_MDS_coord{i}(j,1),mean_MDS_coord{i}(j,2),'ko','MarkerSize',6,'MarkerFaceColor',col);
        text(mean_MDS_coord{i}(j,1)+0.013,mean_MDS_coord{i}(j,2),num2cellstr(ceil(j/2)),'FontSize',14)
    end
end

ind_triu = find(tril(ones(30,30),-1));
model = ones(30,30);
model = model - eye(30,30);
model(16:30,16:30)=0;

ind = [1:2:30,2:2:30];
for i=1:size(coord,1)
    %subplot(2,2,i);
    mat = squareform(mean_D{i});
    mat=mat(ind,ind);
    plot_matrix_simple(mat,sprintf('r = %1.3f (x=%i,y=%i,z=%i)',corr(mat(ind_triu),model(ind_triu),'type','spearman'),coord_mni(i,1),coord_mni(i,2),coord_mni(i,3)));
end


% figure;
% for i=1:size(coord,1)
%     subplot(2,2,i);
%     plot(mean_shape{i}(:,1),mean_shape{i}(:,2),'or','MarkerSize',6,'MarkerFaceColor','r');
%     title(sprintf('%s\n(x=%i,y=%i,z=%i)',roilabel{i}{1},coord_mni(i,1),coord_mni(i,2),coord_mni(i,3)));
%     text(mean_shape{i}(:,1)+0.013,mean_shape{i}(:,2),num2cellstr(1:30),'FontSize',14)
% end