close all;

for i=1:length(full_corrdata)
    
    str = sprintf('Key-frame %i',TIMING_data(i).ID);
    
    y = full_corrdata(i).FIRST_times;
    x = full_corrdata(i).SECOND_times;
    
    mat = full_corrdata(i).corrmat.*(full_corrdata(i).corrmat>0.9);
    if nnz(mat)==0
       error('Tyhjä matriisi!');
    end
    plot_matrix_simple(mat,str,num2cellstr(y),num2cellstr(x));
    xlabel(sprintf('K2 (second showing AKA key-frame), SES %i',full_corrdata(i).FIRST_ses));
    ylabel(sprintf('K1 (initial showing), SES %i',full_corrdata(i).SECOND_ses));
    grid on
    
    hold on
    
    xx = [minind(abs(x-TIMING_data(i).SECOND_time(1))),minind(abs(x-TIMING_data(i).SECOND_time(2)))];
    yy = [minind(abs(y-TIMING_data(i).FIRST_time(1))),minind(abs(y-TIMING_data(i).FIRST_time(2)))];
    
    plot(xx,yy,'go','MarkerSize',10);
    plot(xx,yy,'kx','MarkerSize',10);
    
    saveas(gcf,['KEYFRAME_',num2str(TIMING_data(i).ID),'_matrix.png']);
end