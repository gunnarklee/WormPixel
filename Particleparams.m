%% get Param range for the worm
function Particleparams (resz, Alldata, DateFldrNms, imgfmt, dynamicTH, thresh_hold, PRCbracketSz)
% get particle filter values by asking the user to idenfity the worm in 5 pictures.
WormProps=[];

for W=1:length(DateFldrNms) % for each folder check for filter params
    dirOutputCropPar = dir(fullfile([Alldata, '/',DateFldrNms{W}], 'FltrParams.mat')); %list the images
    
    if size( dirOutputCropPar,1) < 1 % if there are no filter params, then get them
        %>DateFldrNms2 = {dirOutputCropPar.name}'; % get image list
        dirOutputtif = dir(fullfile([Alldata, '/',DateFldrNms{W}], imgfmt)); %list the images
        
        %Get 5 images at equal intervals across iamge set
        QueryImgs=[1:ceil(size(dirOutputtif,1)/6):ceil(size(dirOutputtif,1))];
        for imgnm=QueryImgs %load 5 images to choose particles
            imgtmp=imread([Alldata filesep DateFldrNms{W} filesep dirOutputtif(imgnm).name]); %load the first image in the folder
            imgtmp=imresize(imgtmp, resz);
            if imgnm == 1; BgndMsk=ones(size(imgtmp));end %make empty mask
            if size(imgtmp,3)==3; imgtmp=rgb2gray(imgtmp); end
            
            if strcmpi(dynamicTH, 'y')
                thresh_hold = graythresh(imgtmp);% dynamically determine threshold
            end
            
            %% Get thresholded image
            % allow threshold adjustment in first image
            if imgnm==1
                TH = 'y';
                while strcmpi (TH, 'n') == 0
                    imgBW=imcomplement(im2bw(imgtmp,thresh_hold));% apply threshold
                    figure; imagesc(imgBW);
                    display(thresh_hold)
                    TH=input('rethreshold? - value or "n"','s');
                    if strcmpi(TH, 'n')==0
                        thresh_hold=str2num(TH);%=input('new value')
                        dynamicTH = 'n';
                    end
                end
                display(thresh_hold)
            else
                imgBW=imcomplement(im2bw(imgtmp,thresh_hold));% apply threshold
            end
            
            
            %buld a mask that collects image artifacts
            %BgndMsk=BgndMsk(:,:,1).*imgBW;
            %figure;imshow(BgndMsk);            
                 
            %find particles and get partilce properties
            [imgBWL, F, Image_PropertiesAll] = GetImgProps (imgBW, 'y');
            if imgnm==1; StartPos=[ceil(size(imgBWL(:,:,1), 2)*.5), ceil(size(imgBWL (:,:,1), 1)*.5)]; end
            
            %identify the worm
            mssg='Double click on the worm centorid'
            [pos] = GetPoint(imgBWL, StartPos, mssg);
            StartPos=pos;
            if imgnm==1; FltrParams.StartPos=StartPos; end
            %get the particle with the closest centroid to the selected location
            [row, mindiff]=CloseCentr(pos, Image_PropertiesAll);
            
            %capture list of worm props
            Currhit=Image_PropertiesAll(row,:);
            WormProps=[WormProps;Currhit];
        end
        
        %build and save new parameter spec sheet
        %% BUILD THE NEW PARTICLE FILTERS
        [upper, lower, stats]= GetParamLimits(WormProps(:,3), PRCbracketSz);
        FltrParams.ParticleFilt.LowLim=lower;
        FltrParams.ParticleFilt.UpLim=upper;  %col10 <<NEEDED TO RAise AREA TO <85 for clump; <5 for smallest only
        
        %MAJOR and MINOR AXIS NEED TO BE MORE PERMISSIVE
        [upper, lower, stats]= GetParamLimits(WormProps(:,8), .8);
        FltrParams.ParticleFilt.MajAxL=lower;
        FltrParams.ParticleFilt.MajAxU=upper; %col 8   <<NEEDED TO RAise AREA TO >24 for clump;  <3.3 for smallest only
        
        [upper, lower, stats]= GetParamLimits(WormProps(:,9), .8);
        FltrParams.ParticleFilt.MinAxL=lower;
        FltrParams.ParticleFilt.MinAxU=upper; %col 9   <<NEEDED TO RAise AREA TO >12 for clump;
        
        FltrParams.TotAxL=FltrParams.ParticleFilt.MinAxL+FltrParams.ParticleFilt.MajAxL;
        FltrParams.TotAxU=FltrParams.ParticleFilt.MinAxU+FltrParams.ParticleFilt.MajAxU;
        FltrParams.threshold=thresh_hold
        %FltrParams.BgndMsk=BgndMsk
        
        %% filters
        FiltApp='SZ_Ax_BB';%'all'
        save ([Alldata filesep DateFldrNms{W} filesep 'FltrParams.mat'], 'FltrParams', 'pos'); %save parm -posctr s ellipse
        close all %clean up
    end
end %folder loop
end