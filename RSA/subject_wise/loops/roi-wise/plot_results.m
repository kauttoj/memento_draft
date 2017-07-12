clc;
clear all;
close all;

ROIs_of_interest{1} = [... % 0sec
    120
    74
    176
    183
    179
    200
    170];

ROIs_of_interest{2} = [... % -5sec
    127
    136
    102
    155
    131
    52
    78
    135];

load volume_selection_data.mat

% for subj_ind = 1:length(S);    
%     subj = S{subj_ind};
%     load([subj,'_correlation_results.mat']);    
%    
%     figure('name',sprintf('Subject %s',subj));
%     for i=1:length(ROIs_of_interest)        
%         ind = ROIs_of_interest{i}(:)';        
%         subplot(length(ROIs_of_interest),1,i);hold on;
%         for j=ind
%             plot(DELAYS,all_correlations(:,j));
%         end
%         axis tight;
%         title(sprintf('ROI set %i',i));
%         xlabel('Delay')
%         ylabel('Correlation')
%         box on;
%     end
%     
% end

load('combined_nullvals.mat');

for i=1:length(ROIs_of_interest)
    
    ind = ROIs_of_interest{i}(:)';
    subplot(length(ROIs_of_interest),1,i);hold on;
    
    mean_correlations = zeros(length(DELAYS),length(ind));
    
    for subj_ind = 1:length(S);        
        subj = S{subj_ind};
        load([subj,'_correlation_results.mat']);                
        for j=1:length(ind)
            mean_correlations(:,j) = mean_correlations(:,j) + all_correlations(:,ind(j));
        end        
    end    
    mean_correlations=mean_correlations/length(S);
    
    plot(DELAYS,mean_correlations);
    
    axis tight;
    title(sprintf('ROI set %i',i));
    xlabel('Delay')
    ylabel('Correlation')
    box on;
    
end