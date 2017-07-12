function print_frames(corrdata,ORIGINAL_data,KRON_data,i,j,str,WIDTH,HEIGHT,ses1,ses2)
%PRINT_FRAMES Summary of this function goes here
%   Detailed explanation goes here



fprintf('Best correlation is %.4f is at %.2fs and %.2fs\n', corrdata(i,j),ORIGINAL_data(i).timepoint,KRON_data(j).timepoint);

scrsz = get(groot,'ScreenSize');
handle = figure('Position',[1,scrsz(4),scrsz(3)+10,scrsz(4)/2.2]);

subplot(1,2,1);
image(uint8(ORIGINAL_data(i).cdata));
set(gca,'units','pixels','Position',[0,0,scrsz(3)/2 - 5,HEIGHT]);
box off
set(gca,'XTick',[],'YTick',[]);
aa=title(sprintf('ORIGINAL ses %i: %0.2fs (%s)',ses1,ORIGINAL_data(i).timepoint,sec2min(ORIGINAL_data(i).timepoint)));

subplot(1,2,2);
image(uint8(KRON_data(j).cdata));
set(gca,'units','pixels','Position',[scrsz(3)/2 + 5,0,scrsz(3)/2,HEIGHT]);
box off
set(gca,'XTick',[],'YTick',[]);
bb=title(sprintf('KRON ses %i: %0.2fs (%s)',ses2,KRON_data(j).timepoint,sec2min(KRON_data(j).timepoint)));

annotation('textbox',[0.4,0.97,0.2,eps],'String',sprintf('correlation %.3f',corrdata(i,j)),'LineStyle','none','HorizontalAlig','center');

set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf,[str,'.png']);
%saveas(gcf,[str,'.tif']);

close

end