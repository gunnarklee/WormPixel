%% get Param range for the worm
function MeasureScale (resz, Alldata)
% get particle filter values by asking the user to idenfity the worm in 5 pictures.
dirOutputCropPar = dir(fullfile([Alldata  filesep 'Measure.mat'])); %list the images
if size( dirOutputCropPar,1) < 1 % if there are no filter params, then get them

    calib= input ('Is there a claibration image? (y/n)', 's');
    if strcmpi(calib, 'y');
        [ScaleImage, path, ~]=uigetfile({'*.jpg;*.tif;*.png;*.gif;*.*','All Image Files'}, 'select the scale Image stretch the line to 1mM', Alldata)%InputDir
        ScaleIM=imread([path filesep ScaleImage]);
        ScaleIM=imresize(ScaleIM, resz);
        figure; imshow(ScaleIM);
        h = imline;
        position = wait(h);
        SCALETYPE='(mM)';
        ImgSz=size(ScaleIM)
        %pdist2([position(1,1),position(1,2)],[position(2,1),position(2,2)])%SAME
        pixelPERmm=sqrt((position(1,1)-position(2,1))^2 + (position(1,2)-position(2,2))^2 );
        SCALETYPE='mm'
    else 
        pixelPERmm=1
        ImgSz=[2500,4500]
         SCALETYPE='pixels'
    end
        
    save ([Alldata filesep 'Measure.mat'], 'pixelPERmm',  'ImgSz', 'SCALETYPE'); %save parm -posctr s ellipse
    close all %clean up
end

end