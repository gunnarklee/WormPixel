%function FltrandDisp
%G.kleemann 2/6/08
%version 3 adds optional image supression

% size filter to exclude unwanted vaules (1 or 9..46 and above) -(Image_PropertiesAll(:,3)
% this will be done twice, once for "long" or worm like objects, once for "round" or egg like objects

%% Apply any of the filters and display the results...Integrate a quesiton with a switch statement for testing
%you can combine filters by multipling them

countgood = 'n';
FiltersTried=1;

while strcmpi('n', countgood)
    
    switch FiltApp
        case 'all'
            filter= [BndBxFlt4L.* BndBxFlt4U.* szFltU.* szFltL.*EccFltU.* EccFltL.*MajAxFltU.* MajAxFltL.*MinAxFltU.* MinAxFltL.*ExtFltL.*ExtFltU];
        case 'SZ_Ax_BB'
            filter= [szFltU.* szFltL.*MajAxFltU.* MajAxFltL.*MinAxFltU.* MinAxFltL.*TotAxFltU.*TotAxFltL];    %BndBxFlt4L.* BndBxFlt4U.*
        case 'uin'
            filter= filter; %specified as filter by callng program
    end
    
    clear ('filtVal')
    filtVal=find(filter); %finds row addresses of values

    Img_Propfilt=Image_PropertiesAll(filtVal,:); %new matrix with 0s filtered out
       

%% GET HEAD POSITON STEP
         %HeadPosnLs=[]...$get this from the skeleton

    %% PARTICLE COUNT CHECK
    %IF MORE PARTICLES ARE DETECTED THAN ARE KNOWN TO BE ON PLATE
    %RE-filter with the more stringent set of parameters
    %expect at max 2 times max number

    %% PARTICLE THRESHOLD AND FILTER SWITCHES
    if and(~size(Img_Propfilt, 1) == MinWorms*MaxWormFactor, FiltersTried < MaxFilt);
            %just try one other filter
        countgood = 'n';
        StringentFilter = 'y';
        FiltersTried=FiltersTried+1;
       % OneWorm_CHR7Params  % ** dont overwrite adapted params **reload params for next filtering
    elseif and (~size(Img_Propfilt, 1) == MinWorms*MaxWormFactor, FiltersTried == MaxFilt);
        % if you have tried all the filters and htere are still too many
        % worms, discard the time point
        countgood = 'n'; %skip image display and mark for continue when out of loop
        Skiplog={Skiplog; imageName};
        
        break % next iteration without saving in data
    else% WORMS BELOW COUNT THRESHOLD size(Img_Propfilt, 1)<=MaxWorms*MaxWormFactor
        countgood = 'y';
    end
end

if strcmpi('y', countgood)
    %% COLLECT DATA
    centroidLs= [centroidLs; Img_Propfilt(1,1:2)];% add a centroid
    centrSmthY=smooth(centroidLs, SmoothInt, SmoothMeth); % smooth in one set
    centrSmthY=reshape(centrSmthY, size( centroidLs, 1), size( centroidLs, 2));
        
        if length(centrSmthY) == size(centroidLs, 1)-1; 
        centrSmthY=[centroidLs(1,:);centrSmthY]; % need to add a spacer.
        end
        
    CentrComp= [centroidLs, centrSmthY]; % rebuilding CentrComp each time
    %headposn
    
    Imagesfilt={};
    for ChosenIMAGE=1:length(filtVal);
        Imagesfilt = [Imagesfilt; F(filtVal(ChosenIMAGE,1)).Image];
    end    
    
    %% EXTRACT THE object OF INTEREST
    boundingBox=Img_Propfilt (1, 11:14);
    boundingBox(3)= boundingBox(3)-1;
    boundingBox(4)= boundingBox(4)-1;
    ImgCellPic=double(imcrop(img1, boundingBox));
    ImgCellmsk=double(Imagesfilt{1,1});
    ImgCell=(ImgCellmsk.*ImgCellPic);
    
    
%% extract partilce params
    areaboundingBox= boundingBox(3)* boundingBox(4);
    ImgCellNoZero=ImgCell(ImgCell~=0);
    areacell = (Img_Propfilt(1,3));

        %major/minor axis
        if boundingBox(3)> boundingBox(4) ;
        MajvsMin=boundingBox(3)/ boundingBox(4);
        else
        MajvsMin=boundingBox(4)/ boundingBox(3);
        end
        BBratio=[BBratio; MajvsMin, imageCt];

    numObj=length (Img_Propfilt (:,1));
    
    %% results
    
    if strcmpi(allow_img, 'y')
        %HISTOGRAM filtered values
        figure;
        title ('all filters applied')
        subplot (3,4,1); hist (Img_Propfilt  (:,5)); title ('Eccentr 1=straight')
        subplot (3,4,2); hist (Img_Propfilt  (:,8)); title ('MajAxis')
        subplot (3,4,3); hist (Img_Propfilt (:,9)); title ('MinAxis')
        subplot (3,4,4); hist (Img_Propfilt (:,10)); title ('Area')
        subplot (3,4,5); hist (Img_Propfilt  (:,15)); title ('Extent')
        subplot (3,4,6); hist (Img_Propfilt (:,16)); title ('equivalentDiameter')
        subplot (3,4,7); hist (Img_Propfilt  (:,14)); title ('EulerNumber')
        subplot (3,4,8); hist (Img_Propfilt (:,18)); title ('Perimiter/area')
        subplot (3,4,9); hist (Img_Propfilt  (:,20)); title ('maj_min_extnt')
        subplot (3,4,10); hist (Img_Propfilt  (:,14)); title ('BOundingBox4(col-14')
        
        % DISPLAY RESULTS
        close all
        figure; imshow(imgBW); %axis equal;
        hold on;
        plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'w*')
        
        figure ('position', scrsz); imagesc(img1); %axis equal;
        hold on;
        plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*')
        
        %ADD PARTICLE NUMBER
        for rowN = 1:length(Img_Propfilt(:,1))
            text(Img_Propfilt(rowN,6),...
                Img_Propfilt(rowN,7),...
                num2str(Img_Propfilt(rowN,4)),'Color', 'r', 'FontSize', 12)
        end
        
        %
        figure ('position', scrsz); imagesc(img-2); %axis equal;
        hold on;
        plot (Image_PropertiesAll(:,6), Image_PropertiesAll(:,7),'wO')
        %plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'b*')
        %ADD PARTICLE NUMBER
    end
    
    
    Image_PropertiesList = Image_PropertiesAll; ctX = 6; ctY = 7; TxtCol = 4; NumCol = 'r'; FntSz = 6;
    if strcmpi ('y', allow_img);
        %cd(CodeDir)
        PrintObjParamsAug11
    end
    
    %cd([Alldata '/' TrialName 'RUNfinal'])
    %% FILTERED FIGURE FOR SCORING DIAGNOSIS
    
    if strcmpi(ProofingImgs, 'y')
        nameProof= [DateFldrNms{W},'MultView',imageName,'.pdf'];
        ProofImages(scrsz, img1, Img_Propfilt, img, nameProof) %make function
    end %end proofing images loop
      
    %% CAPTURE AS STACKED TIFF
    if strcmpi(DataCapMode, 'StackGiff')
        
          % Append to stack of images
       %Filename = [RUNfiltDir, filesep, DateFileNms{Pair}(1:end-9), 'stack.gif'];
        Filename = [RUNfinalDir filesep DateFldrNms{W},'stack.gif'];
        measure = ['BB-ratio',num2str(MajvsMin)];%
        textls={'Thrashing in Buffer'; imageName; [num2str(timeintv), 'secs']; measure};
        Nm2=strrep(Filename, '_', '-');
        saveImageToStack(uint8(img1), Filename, ...
            'title', 'Image', ...
            'image_name', Nm2, ...
            'scale_image', false, ...
            'display_image', 'on',...
            'CentrComp', CentrComp,...
            'boundingBox', boundingBox,...
            'plotcol',[1,2,3,4]); %proofingImgVIS
     end

end

close all


 