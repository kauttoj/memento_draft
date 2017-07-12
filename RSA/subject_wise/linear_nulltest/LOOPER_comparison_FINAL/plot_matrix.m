function handle = plot_matrix(mat,tit,labels_x,labels_y)

M = 4*140;

N_colors = 2*M + 1;
%center = ceil(N_colors/2);

if iscell(mat)
   mat_plot=mat{2};
   mat=mat{1};    
else
   mat_plot=mat;
end

mat = double(mat);

Lx=size(mat,1);
Ly=size(mat,2);

val = mat(:);

ma = max(val);
mi =  min(val);

if sign(ma) == -sign(mi)
    colormapp1 = autumn(M);
    colormapp2 = flipud(winter(M));
    colormapp = [colormapp1;1,1,1;colormapp2];
    
    a=linspace(ma,0,6);
    b=linspace(0,mi,6);
    a = [a(1:5),0,b(2:6)];
    for i=1:length(a)
        colorbar_labels{i} = sprintf('%1.2f',a(i));
    end
    
else
    colormapp = jet(N_colors-1);
    colormapp = flipud(colormapp);
    colormapp(end+1,:)=[1,1,1];
    
    a=linspace(ma,0,11);
    for i=1:length(a)
        colorbar_labels{i} = sprintf('%1.0f',a(i));
    end
end

orig_mat = mat_plot;

ind = orig_mat<=0;
if mi<0
    mat_plot(ind) = round(M*((orig_mat(ind)/mi)))+1+M;
else
    mat_plot(ind) = round(M*((orig_mat(ind))))+1+M;
end

ind = orig_mat>0;
mat_plot(ind) = round(M*(1-(orig_mat(ind)/ma)))+1;

v = mat_plot(ind);
v(v==M+1)=M;
mat_plot(ind)=v;

%mat(mat_zeros)=0;
% colorbar_min = mi;
% colorbar_max = ma;

colorbar_ticks = round(linspace(1,size(colormapp,1),11));


totalmap = colormapp;

if length(labels_y)>70
    handle = figure('Position',[75   1   906   1200],'Name',tit);
else
    handle = figure('Position',[75   387   906   850],'Name',tit);
end

orig_axis = axes('Position',[0.5 0.06 0.37 0.93],'Visible','on');

set(handle, 'PaperPositionMode', 'auto');

annotation('textbox',[0.1 0.94 0.8 0.05],'string',tit,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','none')

imagesc(1:size(mat_plot,2),1:size(mat_plot,1),mat_plot);

%orig_axis = gca;

ax1 = axes('Position',[0 0 1 1],'Visible','off');

if nargin>=6 

  un = unique(super_label_id);
  ma = length(super_label_id);
  for id = un
      k=find(super_label_id==id);
      
      set(handle,'CurrentAxes',orig_axis)
      [xaf,yaf] = ds2nfu([0,Lx],ma-[k(1),k(1)]+1.5);
      annotation('line',[xaf(1)*0.8,xaf(2)*1.01],[yaf(1),yaf(2)]);
      [xaf,yaf] = ds2nfu([0,Lx],ma-[k(end),k(end)]+0.5);
      annotation('line',[xaf(1)*0.8,xaf(2)*1.01],[yaf(1),yaf(2)]);
      
      [xaf,yaf] = ds2nfu(0,ma-(k(1)+k(end))/2 + 0.5 + 1);      
      set(handle,'CurrentAxes',ax1)
      text(xaf*0.3,yaf,super_label{id},'FontSize',12,'Rotation',0);%,'HorizontalAlignment','left');
      %annotation('textbox',[xaf*0.3,yaf,0.05,0.05],'String',{super_label{id}},'FontSize',12, 'FitBoxToText','on','EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle')     
      
      set(handle,'CurrentAxes',orig_axis)
      [xaf,yaf] = ds2nfu([k(1),k(1)]-0.5,[0,Ly]);
      annotation('line',[xaf(1),xaf(2)],[yaf(1)*0.8,yaf(2)*1.01]);
      [xaf,yaf] = ds2nfu([k(end),k(end)]+0.5,[0,Ly]);
      annotation('line',[xaf(1),xaf(2)],[yaf(1)*0.8,yaf(2)*1.01]);      
      
      [xaf,yaf] = ds2nfu((k(1)+k(end))/2+0.5 - 0,0);
      set(handle,'CurrentAxes',ax1)
      text(xaf,yaf*0.1,super_label{id},'FontSize',12,'Rotation',90);%,'HorizontalAlignment','left');
      
%      annotation('textbox',[xaf,yaf*0.3,0.05,0.05],'String',{super_label{id}},'FontSize',10, 'FitBoxToText','on','EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle')     
      
  end
end
set(handle,'CurrentAxes',orig_axis)
if nargin>3 && ~isempty(labels_x) && ~isempty(labels_y)
            
    set(gca,'YTick',1:length(labels_y));
    
    %set(gca,'YTickLabel',labels);
    
    set(gca,'XTick',1:length(labels_x));
    
    %[hx,hy]=format_ticks(gca,labels_x,labels_y,[],[],0,0);
    %labels_x = shorten_labels(labels_x);
    %labels_y = shorten_labels(labels_y);
    
    labels_x = shorten_labels(labels_x);
    
%    set(gca,'XTickLabel',labels_x,'fontsize',16)

    %xticklabel_rotate(1:size(mat,2),0,labels_x,'fontsize',16);
        
    if length(labels_y)>70
        set(gca,'YTickLabel',labels_y,'fontsize',8)
            xticklabel_rotate(1:length(labels_x),90,labels_x,'fontsize',10);

    else
        set(gca,'YTickLabel',labels_y,'fontsize',14)
            xticklabel_rotate(1:length(labels_x),90,labels_x,'fontsize',14);

    end
    
%     set(hy,'fontsize',12)
%     set(hx,'fontsize',12)

elseif  nargin>2 && ~isempty(labels_y)
    set(gca,'YTick',1:length(labels_y));
    %set(gca,'YTickLabel',labels);
    
    set(gca,'XTick',1);
    %xticklabel_rotate(1:size(mat,2),90,labels);
    labels_x = shorten_labels(labels_x);
    
    set(gca,'XTickLabel',labels_x,'fontsize',12)
    %set(gca,'YTickLabel',labels_y,'fontsize',12)
    
    %set(hy,'fontsize',8)
%set(hx,'fontsize',8)
else
    set(gca,'YTick',[]);    
    set(gca,'XTick',[]);      
end

set(gca,'YDir','Normal')

%axis square;
colormap(totalmap);
%imshow(mat,totalmap, 'InitialMagnification','fit');

%set(gca,'YTick',[]);
%set(gca,'XTick',[]);
%box on;
set(gca,'YDir','Reverse')
%set(gca,'YGrid','on')

%xlabel('Raw','FontSize',16)
%ylabel('Convolved & downsampled','FontSize',16)

%hold on;




%ax = axes('Position', [0.05 0.3 0.9 0.4], 'Visible', 'off');
h = colorbar('FontSize',12,'Position',[0.88,0.4,0.03,0.23]);
lim = get(h,'ylim');
%set(h, 'ylim',[1,N_colors]);
scal = (max(lim)-min(lim))/(max(colorbar_ticks)-min(colorbar_ticks));
set(h, 'YTick',min(colorbar_ticks)+(colorbar_ticks-min(colorbar_ticks))*scal);
set(h,'YTickLabel',(colorbar_labels));
set(h,'YDir','reverse')


end


% scale_subplots(1.12,[-0.05,-0.05])
%
% figure;
% plot(sort(val));
% xlabel('value number');
% ylabel('Z-score');
% axis tight;
function new_labs = shorten_labels(labs)

    
    if length(labs)>50
        inds = round(linspace(1,length(labs),10));
        for i=1:length(labs)
                        
            if ismember(i,inds)
                new_labs{i}=num2str(labs(i));
            else
                new_labs{i}='';
            end
            
        end
    else
        new_labs = labs;
    end

    new_labs = num2cellstr(new_labs);    

end

%%% oma koodi loppuu tähän

function hText = xticklabel_rotate(XTick,rot,varargin)
%hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
%
% Syntax: xticklabel_rotate
%
% Input:    
% {opt}     XTick       - vector array of XTick positions & values (numeric) 
%                           uses current XTick values or XTickLabel cell array by
%                           default (if empty) 
% {opt}     rot         - angle of rotation in degrees, 90� by default
% {opt}     XTickLabel  - cell array of label strings
% {opt}     [var]       - "Property-value" pairs passed to text generator
%                           ex: 'interpreter','none'
%                               'Color','m','Fontweight','bold'
%
% Output:   hText       - handle vector to text labels
%
% Example 1:  Rotate existing XTickLabels at their current position by 90�
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45� and change
% font size
%    xticklabel_rotate([],45,[],'Fontsize',14)
%
% Example 3:  Set the positions of the XTicks and rotate them 90�
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45� without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90� at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Note : you can not RE-RUN xticklabel_rotate on the same graph. 
%



% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%                       Arbitrary angle rotation
%                       Output of text handles
%                       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size 
%                           and keep text on plot
%                           (handles small window resizing after, but not well due to proportional placement with 
%                           fixed font size. To fix this would require a serious resize function)
%                       Uses current XTick by default
%                       Uses current XTickLabel is different from XTick values (meaning has been already defined)

% Brian FG Katz
% bfgkatz@hotmail.com
% 23-05-03
% Modified 03-11-06 after user comment
%	Allow for exisiting XTickLabel cell array
% Modified 03-03-2006 
%   Allow for labels top located (after user comment)
%   Allow case for single XTickLabelName (after user comment)
%   Reduced the degree of resizing
% Modified 11-jun-2010
%   Response to numerous suggestions on MatlabCentral to improve certain
%   errors.

% Other m-files required: cell2mat
% Subfunctions: none
% MAT-files required: none
%
% See also: xticklabel_rotate90, TEXT,  SET

% Based on xticklabel_rotate90
%   Author: Denis Gilbert, Ph.D., physical oceanography
%   Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%   email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%   February 1998; Last revision: 24-Mar-2003

% check to see if xticklabel_rotate has already been here (no other reason for this to happen)
if isempty(get(gca,'XTickLabel')),
    error('xticklabel_rotate : can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased')  ;
end

% if no XTickLabel AND no XTick are defined use the current XTickLabel
%if nargin < 3 & (~exist('XTick') | isempty(XTick)),
% Modified with forum comment by "Nathan Pust" allow the current text labels to be used and property value pairs to be changed for those labels
if (nargin < 3 || isempty(varargin{1})) & (~exist('XTick') | isempty(XTick)),
	xTickLabels = get(gca,'XTickLabel')  ; % use current XTickLabel
	if ~iscell(xTickLabels)
		% remove trailing spaces if exist (typical with auto generated XTickLabel)
		temp1 = num2cell(xTickLabels,2)         ;
		for loop = 1:length(temp1),
			temp1{loop} = deblank(temp1{loop})  ;
		end
		xTickLabels = temp1                     ;
	end
varargin = varargin(2:length(varargin));	
end

% if no XTick is defined use the current XTick
if (~exist('XTick') | isempty(XTick)),
    XTick = get(gca,'XTick')        ; % use current XTick 
end

%Make XTick a column vector
XTick = XTick(:);

if ~exist('xTickLabels'),
	% Define the xtickLabels 
	% If XtickLabel is passed as a cell array then use the text
	if (length(varargin)>0) & (iscell(varargin{1})),
        xTickLabels = varargin{1};
        varargin = varargin(2:length(varargin));
	else
        xTickLabels = num2str(XTick);
	end
end    

if length(XTick) ~= length(xTickLabels),
    error('xticklabel_rotate : must have same number of elements in "XTick" and "XTickLabel"')  ;
end

%Set the Xtick locations and set XTicklabel to an empty string
set(gca,'XTick',XTick,'XTickLabel','')

if nargin < 2,
    rot = 90 ;
end

% Determine the location of the labels based on the position
% of the xlabel
hxLabel = get(gca,'XLabel');  % Handle to xlabel
xLabelString = get(hxLabel,'String');

% if ~isempty(xLabelString)
%    warning('You may need to manually reset the XLABEL vertical position')
% end

set(hxLabel,'Units','data');
xLabelPosition = get(hxLabel,'Position');
y = xLabelPosition(2);

%CODE below was modified following suggestions from Urs Schwarz
y=repmat(y,size(XTick,1),1);
% retrieve current axis' fontsize
fs = get(gca,'fontsize');

% Place the new xTickLabels by creating TEXT objects
hText = text(XTick, y, xTickLabels,'fontsize',fs);

% Rotate the text objects by ROT degrees
%set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
% Modified with modified forum comment by "Korey Y" to deal with labels at top
% Further edits added for axis position
xAxisLocation = get(gca, 'XAxisLocation');  
if strcmp(xAxisLocation,'bottom')  
    set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})  
else  
    set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:})  
end

% Adjust the size of the axis to accomodate for longest label (like if they are text ones)
% This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
% This approach keeps the right side of the graph at the same place 

set(get(gca,'xlabel'),'units','data')           ;
    labxorigpos_data = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'ylabel'),'units','data')           ;
    labyorigpos_data = get(get(gca,'ylabel'),'position')  ;
set(get(gca,'title'),'units','data')           ;
    labtorigpos_data = get(get(gca,'title'),'position')  ;

set(gca,'units','pixel')                        ;
set(hText,'units','pixel')                      ;
set(get(gca,'xlabel'),'units','pixel')          ;
set(get(gca,'ylabel'),'units','pixel')          ;

origpos = get(gca,'position')                   ;

% textsizes = cell2mat(get(hText,'extent'))       ;
% Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given. 
x = get( hText, 'extent' );  
if iscell( x ) == true  
    textsizes = cell2mat( x ) ;  
else  
    textsizes = x;  
end  

largest =  max(textsizes(:,3))                  ;
longest =  max(textsizes(:,4))                  ;

laborigext = get(get(gca,'xlabel'),'extent')    ;
laborigpos = get(get(gca,'xlabel'),'position')  ;

labyorigext = get(get(gca,'ylabel'),'extent')   ;
labyorigpos = get(get(gca,'ylabel'),'position') ;
leftlabdist = labyorigpos(1) + labyorigext(1)   ;

% assume first entry is the farthest left
leftpos = get(hText(1),'position')              ;
leftext = get(hText(1),'extent')                ;
leftdist = leftpos(1) + leftext(1)              ;
if leftdist > 0,    leftdist = 0 ; end          % only correct for off screen problems

% botdist = origpos(2) + laborigpos(2)            ;
% newpos = [origpos(1)-leftdist longest+botdist origpos(3)+leftdist origpos(4)-longest+origpos(2)-botdist]  
%
% Modified to allow for top axis labels and to minimize axis resizing
if strcmp(xAxisLocation,'bottom')  
    newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
            origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
            origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
            origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))];
else
    newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
            origpos(2) ...
            origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
            origpos(4)-(longest)+get(gca,'FontSize')];
end
set(gca,'position',newpos);

% readjust position of text labels after resize of plot
set(hText,'units','data');
for loop= 1:length(hText),
    set(hText(loop),'position',[XTick(loop), y(loop)])  ;
end

% adjust position of xlabel and ylabel
laborigpos = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0])   ;

% switch to data coord and fix it all
set(get(gca,'ylabel'),'units','data')                   ;
set(get(gca,'ylabel'),'position',labyorigpos_data)      ;
set(get(gca,'title'),'position',labtorigpos_data)       ;

set(get(gca,'xlabel'),'units','data')                   ;
    labxorigpos_data_new = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)])   ;


% Reset all units to normalized to allow future resizing
set(get(gca,'xlabel'),'units','normalized')          ;
set(get(gca,'ylabel'),'units','normalized')          ;
set(get(gca,'title'),'units','normalized')          ;
set(hText,'units','normalized')                      ;
set(gca,'units','normalized')                        ;

if nargout < 1,
    clear hText
end

end

function out = imoverlay(in, mask, color)
%IMOVERLAY Create a mask-based image overlay.
%   OUT = IMOVERLAY(IN, MASK, COLOR) takes an input image, IN, and a binary
%   image, MASK, and produces an output image whose pixels in the MASK
%   locations have the specified COLOR.
%
%   IN should be a grayscale or an RGB image of class uint8, uint16, int16,
%   logical, double, or single.  If IN is double or single, it should be in
%   the range [0, 1].  If it is not in that range, you might want to use
%   mat2gray to scale it into that range.
%
%   MASK should be a two-dimensional logical matrix.
%
%   COLOR should be a 1-by-3 vector of values in the range [0, 1].  [0 0 0]
%   is black, and [1 1 1] is white.
%
%   OUT is a uint8 RGB image.
%
%   Examples
%   --------
%   Overlay edge detection result in green over the original image.
%       
%       I = imread('cameraman.tif');
%       bw = edge(I, 'canny');
%       rgb = imoverlay(I, bw, [0 1 0]);
%       imshow(rgb)
%
%   Treating the output of peaks as an image, overlay the values greater than
%   7 in red.  The output of peaks is not in the usual grayscale image range
%   of [0, 1], so use mat2gray to scale it.
%
%       I = peaks;
%       mask = I > 7;
%       rgb = imoverlay(mat2gray(I), mask, [1 0 0]);
%       imshow(rgb, 'InitialMagnification', 'fit')

%   Steven L. Eddins
%   Copyright 2006-2012 The MathWorks, Inc.

% If the user doesn't specify the color, use white.
DEFAULT_COLOR = [1 1 1];
if nargin < 3
    color = DEFAULT_COLOR;
end

% Force the 2nd input to be logical.
mask = (mask ~= 0);

% Make the uint8 the working data class.  The output is also uint8.
in_uint8 = im2uint8(in);
color_uint8 = im2uint8(color);

% Initialize the red, green, and blue output channels.
if ndims(in_uint8) == 2
    % Input is grayscale.  Initialize all output channels the same.
    out_red   = in_uint8;
    out_green = in_uint8;
    out_blue  = in_uint8;
else
    % Input is RGB truecolor.
    out_red   = in_uint8(:,:,1);
    out_green = in_uint8(:,:,2);
    out_blue  = in_uint8(:,:,3);
end

% Replace output channel values in the mask locations with the appropriate
% color value.
out_red(mask)   = color_uint8(1);
out_green(mask) = color_uint8(2);
out_blue(mask)  = color_uint8(3);

% Form an RGB truecolor image by concatenating the channel matrices along
% the third dimension.
out = cat(3, out_red, out_green, out_blue);
end

function new_img = upsample(img,factor)

if factor == 1
    new_img = img;
else
    rows = size(img,1);
    cols = size(img,2);
    layers = size(img,3);
    if layers > 1
        new_img = zeros(rows*factor,cols*factor,layers);
        a = (1:factor:(rows*factor - 1));
        b = (1:factor:(cols*factor - 1));
        for l = 1:layers
            new_img(a,b,l)=img;
            for k=1:(factor-1),
                new_img(a,b,l)=img;
                new_img(a + k,b + k,l)=img;
                new_img(a,b + k,l)=img;
                new_img(a+k,b,l)=img;
            end
        end
    else
        new_img = zeros(rows*factor,cols*factor);
        a = (1:factor:(rows*factor - 1));
        b = (1:factor:(cols*factor - 1));
        new_img(a,b)=img;
        for k=1:(factor-1),
            new_img(a,b)=img;
            new_img(a + k,b + k)=img;
            new_img(a,b + k)=img;
            new_img(a+k,b)=img;
        end
    end
end

end