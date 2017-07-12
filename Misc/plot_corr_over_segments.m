clear m
clear IDs
close all;

for i=1:length(full_corrdata)
    
    try
    
    ID = full_corrdata(i).ID;
    IDs(i) = ID;
    plot_timecourses(ID,full_corrdata(i).total_corr,full_corrdata(i).ORIGINAL_time,full_corrdata(i).KRON_time);
    set(gca,'YLim',[0.5,1]);
    a=full_corrdata(i).total_corr;
    a(a==0)=[];
    a(isnan(a))=[];
    m(i)=mean(a);
    %saveas(gcf,sprintf('ID_%i_correlation_ts.png',ID));
    pause(3);
    catch err
        
    end    
        
end

bar(m)
set(gca,'XTick',1:length(m));
set(gca,'XTickLabel',num2cellstr(IDs));
axis tight;
xlabel('Sequence ID');
ylabel('Mean correlation');
    