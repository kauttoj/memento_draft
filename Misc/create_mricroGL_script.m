function filename = create_mricroGL_script(filename,input_image,savepath)
% USAGE: filename = create_mricroGL_script(filename,input_image,savepath)

if ~iscell(input_image)
    a{1}=input_image;
    input_image=a;
end

if isempty(strfind(filename,filesep))
    filename = [pwd,filesep,filename];
end

if nargin<3 || isempty(savepath) 
    savepath=pwd;
end

if ~exist(savepath,'dir')
    error('outpath does not exist!')
end

if ~strcmp(savepath(end),filesep)
    savepath(end+1)=filesep;
end

fprintf('\n----- 3D slice creator for mricroGL ver. 2012 ------\n');

fileID = fopen([filename,'.gls'],'w');

try
    
    fprintf(fileID,'BEGIN\n');
    fprintf(fileID,'RESETDEFAULTS;\n');
    fprintf(fileID,'CAMERADISTANCE(0.8);\n');
    fprintf(fileID,'LOADIMAGE(''mni152_2009bet'');\n');
    fprintf(fileID,'BACKCOLOR(255, 255, 255);\n');
    fprintf(fileID,'OVERLAYLOADSMOOTH(true);\n');
    fprintf(fileID,'COLORBARVISIBLE(false);\n');
    fprintf(fileID,'EDGEENHANCE(25, 75);\n');
    fprintf(fileID,'OVERLAYTRANSPARENCYONBACKGROUND(10);\n');
    fprintf(fileID,'OVERLAYTRANSPARENCYONOVERLAY(-1);\n');
    
    for i=1:length(input_image)
        fprintf(fileID,'OVERLAYLOAD(''%s'');\n',fullpath(input_image{i}));
    end
    
    fprintf(fileID,'AZIMUTHELEVATION(0,0); // takaa\n');
    fprintf(fileID,'SAVEBMP(''%stakaa.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(0,180); // p??lt?\n');
    fprintf(fileID,'SAVEBMP(''%spaalta.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(180,0); // edest?\n');
    fprintf(fileID,'SAVEBMP(''%sedesta.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(90,0); // vasen\n');
    fprintf(fileID,'SAVEBMP(''%svasen.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(-90,0); // oikea\n');
    fprintf(fileID,'SAVEBMP(''%soikea.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(-140,35); // oik etuviisto\n');
    fprintf(fileID,'SAVEBMP(''%soik_etuviisto.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(140,35); // vas etuviisto\n');
    fprintf(fileID,'SAVEBMP(''%svas_etuviisto.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(40,35); // vas takaviisto\n');
    fprintf(fileID,'SAVEBMP(''%svas_takaviisto.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(-40,35); // oik takaviisto\n');
    fprintf(fileID,'SAVEBMP(''%soik_takaviisto.png'');\n',savepath);
    
    fprintf(fileID,'AZIMUTHELEVATION(180,-90); // alta\n');
    fprintf(fileID,'SAVEBMP(''%salta.png'');\n',savepath);
    
    fprintf(fileID,'CUTOUT(0,0,0,0.5,1,1);\n');
    fprintf(fileID,'AZIMUTHELEVATION(90,0); // oikea sisä\n');
    fprintf(fileID,'SAVEBMP(''%soik_sisa.png'');\n',savepath);
    
    fprintf(fileID,'CUTOUT(0.5,0,0,1,1,1);\n');
    fprintf(fileID,'AZIMUTHELEVATION(-90,0); // vasen sisä\n');
    fprintf(fileID,'SAVEBMP(''%svas_sisa.png'');\n',savepath);
    
    fprintf(fileID,'END.\n');
    
    fclose(fileID);
    
    fprintf('Following script was succesfully created:\n%s\n',filename,'.gls');    
    fprintf('Hint: Use ''image_autocrop'' function to automatically crop images after running the GLS script.\n\n');
    
catch err
    fclose(fileID);
    fprintf('FAILED, reason: %s',err.message);
end

end

function outpath = fullpath(inpath)

outpath = which(inpath);

end

