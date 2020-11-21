function sint_SaveStack_v2(stk,class,path,regname,tifname)
% stk: variable to be saved as a multi-page tiff
% class: uint8 or uint16 depending on numeric input value of 8 or 16
% path: save directory
% regname: name of the new image file to be written
% tifname: name of the original image file, for reading attributes

% note: this function will delete previous instances of the same full file
% path.

if class==8
    stk = uint8(stk);
elseif class==16
    stk = uint16(stk);
end
warning('off','MATLAB:DELETE:FileNotFound');

%%% OLD METHOD: replaced Nov 2017    
%     tic
%     for t0 = 1:size(stk,3)
%         imwrite(stk(:,:,t0),f,'tiff','Compression','none','WriteMode','append');
%     end
% toc
        
%%% NEW METHOD: using direct access to libtiff through Tiff class

% First open the original tiff file as a Tiff object 
im0 = [path tifname];
w = warning('off','all');
y = Tiff(im0,'r+'); % Open object in read mode
warning(w);

% Extract its properties
% (Apart from size of image, which may have changed)
tagstruct.ImageLength=size(stk,1);
tagstruct.ImageWidth=size(stk,2);
tagstruct.Photometric=getTag(y,'Photometric');
tagstruct.SubFileType=getTag(y,'SubFileType');
tagstruct.RowsPerStrip=getTag(y,'RowsPerStrip');
tagstruct.BitsPerSample=getTag(y,'BitsPerSample');
maxstripsize = 8*1024;
tagstruct.RowsPerStrip = ceil(maxstripsize/(size(stk,2)*(tagstruct.BitsPerSample/8)*size(stk,3))); % http://www.awaresystems.be/imaging/tiff/tifftags/rowsperstrip.html
tagstruct.Compression=getTag(y,'Compression');
tagstruct.SampleFormat=getTag(y,'SampleFormat');
tagstruct.SamplesPerPixel= getTag(y,'SamplesPerPixel');
tagstruct.PlanarConfiguration= getTag(y,'PlanarConfiguration');
tagstruct.ImageDescription= getTag(y,'ImageDescription');
tagstruct.Orientation= getTag(y,'Orientation');

close(y); clear y;


% Now create new Tiff object for Registered stack
im1 = fullfile(path,regname); delete(im1);
x = Tiff(im1,'w'); % open in write mode

% Write stack frames to object one by one
for t0 = 1:size(stk,3)    
    setTag(x,tagstruct);      
    x.write(stk(:,:,t0));
    if t0 ~= size(stk,3)
        writeDirectory(x);
    end
end

% Close the object (which saves it to disk)
x.close();  

disp(['Stack ' regname ' saved as a multipage tiff file in ' path '.']);

end