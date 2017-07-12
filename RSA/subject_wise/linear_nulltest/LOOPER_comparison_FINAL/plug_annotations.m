function plug_annotations(handle,table,symbol,fontsize)

a = get(handle,'CurrentAxes');
X = get(a,'XLim');
Y = get(a,'YLim');

siz = size(table);
dx = diff(X)/siz(2);
dy = diff(Y)/siz(1);

for row = 1:size(table,1)
    for col = 1:size(table,2)
        if table(row,col)~=0
            %[xaf,yaf] = ds2nfu(X(1)+col*dx + dx/2,Y(1)+row*dy + dy/2);
            text(col,row,symbol,'FontSize',fontsize,'VerticalAlignment','middle','HorizontalAlignment','center','interpreter','latex');
        end
    end
end

