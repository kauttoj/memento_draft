function image_autocrop(images,background,CUT)
% USAGE: image_autocrop(images,background)
% images = cell of image paths
% background = bacground color to be cropped (w or b)
% CUT = [vasen,ala,oikea,ylä]

if strcmp(images,'all')
   
    i1 = dir('*.png');
    i2 = dir('*.jpg');
    i3 = dir('*.bmp');
    
    images = {};
    
    k=0;
    for i=1:length(i1)
        k=k+1;
        images{k}=i1(i).name;
    end
    for i=1:length(i2)
        k=k+1;
        images{k}=i2(i).name;       
    end
    for i=1:length(i3)
        k=k+1;
        images{k}=i3(i).name;
    end
    
    pois = [];
    for i=1:length(images)        
        if ~isempty(strfind(images{i},'_CROPPED'))
            pois(end+1)=i;
        end
    end
    images(pois)=[];     
    
    fprintf('Found %i images in this folder, cropping them all\n',length(images));
    
end

if ~iscell(images)
    a{1}=images;
    images=a;
end

if nargin<2
    background=255;
else
   if strcmp(background,'w')
       background=255;
   elseif strcmp(background,'b')
       background=1;
   else
       error('Unknown background color (to be cropped), use w or b\n');
   end
end

for i=1:length(images)
    A = imread(images{i});
    
    if nargin==3
        if length(CUT)==4
            fprintf('Cutting edges: left=%i, bottom=%i, right=%i, up=%i\n',CUT(1),CUT(2),CUT(3),CUT(4));
            A=A((1+CUT(4)):(end-CUT(2)),(CUT(1)+1):(end-CUT(3)),:);
        end
    end        
    
    B = (A==background);
    B = sum(B,3);
    B = (B==3);
    
    a=(sum(B)==size(A,1));
    left=find(a==0,1,'first')-1;
    right=find(a==0,1,'last')+1;
    
    a=(sum(B')==size(A,2));
    up=find(a==0,1,'first')-1;
    down=find(a==0,1,'last')+1;
    
    if isempty(left) || left<1
        left=1;
    end
    if isempty(right) || right>size(A,2)
        right=size(A,2);
    end
    if isempty(up) || up<1
        up=1;
    end   
    if isempty(down) || down>size(A,1)
        down=size(A,1);
    end        
    
    A=A(up:down,left:right,:);   
    
    [a,b,c]=fileparts(images{i});
    if isempty(a)
        a=pwd;
    end
    
    imwrite(A,[a,filesep,b,'_CROPPED',c]);
    
end


end

