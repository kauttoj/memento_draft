function plot_timecourses(ID,corrdata,time1,time2 )
%PLOT_TIMECOURSES Summary of this function goes here
%   Detailed explanation goes here

%figure('Position',[443   327   963   498]);

bad = isnan(corrdata);
if nnz(bad)>0
    warning('There are %i NaN values!',nnz(bad));
end
corrdata(bad)=0;

plot(time1,corrdata,'b.-');
%axis tight;
d = (max(time1)-min(time1))/20;
set(gca,'Xlim',[min(time1)-d,max(time1)+d]);
set(gca,'Ylim',[0,1]);
title(sprintf('Memento segment timing check (segment ID %i)',ID));

ind = round(linspace(1,length(time1),10));
val1 = time1(ind);
val2 = time2(ind);
for i=1:length(val1),
    str{i}=sprintf('%s%s%s',sec2min(val1(i),true),'\newline',sec2min(val2(i),true));
end   

ax = gca;

ax.XTick = val1;
ax.XTickLabel =str;

grid on

xlabel 'Time [s]';
ylabel 'Frame correlation (RGB)';

end


