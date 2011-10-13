%function FltrandDisp
%G.kleemann 2/6/08
%version 3 adds optional image supression


% size filter to exclude unwanted vaules (1 or 9..46 and above) -(Image_PropertiesAll(:,3)
% this will be done twice, once for "long" or worm like objects, once for "round" or egg like objects 

%% Apply any of the filters and display the results...Integrate a quesiton with a switch statement for testing 
%you can combine filtes by multipling them
countgood = 'n'
FiltersTried=1
 cd(Alldata)
          %>>>if strcmpi(ProofingImgs, 'y')
          mkdir([TrialName 'RUNfinal'])
          %>mkdir([TrialName 'ScoredPics'])
          mkdir([Alldata '/' TrialName 'SubtrPics'])
          cd([Alldata '/' TrialName 'RUNfinal'])


%%

while strcmpi('n', countgood)

    switch FiltApp 
        case 'no'  
            filter = ones (size (Image_PropertiesAll,1));   
            filter= filter(:,1)    
        case 'all'    
            filter= [BndBxFlt4L.* BndBxFlt4U.* szFltU.* szFltL.*EccFltU.* EccFltL.*MajAxFltU.* MajAxFltL.*MinAxFltU.* MinAxFltL.*ExtFltL.*ExtFltU];
        case 'SZ_Ax_BB'    
            filter= [szFltU.* szFltL.*MajAxFltU.* MajAxFltL.*MinAxFltU.* MinAxFltL];    %BndBxFlt4L.* BndBxFlt4U.*
        case 'sz'
            filter= [szFltU.* szFltL];
        case 'ec'
            filter= [EccFltU.* EccFltL];
        case 'majax'
            filter= [MajAxFltU.* MajAxFltL];
        case 'minax'
           filter= [MinAxFltU.* MinAxFltL];
        case 'ext'
            filter= [ExtFltL.*ExtFltU];
        case 'uin'       
           filter= filter; %specified as filter by callng program
    end

clear ('filtVal')
filtVal=find(filter); %finds row addresses of values
    Img_Propfilt=Image_PropertiesAll(filtVal,:); %new matrix with 0s filtered out

%% PARTICLE COUNT CHECK
%IF MORE PARTICLES ARE DETECTED THAN ARE KNOWN TO BE ON PLATE
%RE-filter with the more stringent set of parameters 
%expect at max 2 times max number

%MaxFilt=2
%MaxWormFactor=2.1

    if and (size(Img_Propfilt, 1)>MaxWorms*MaxWormFactor, FiltersTried < MaxFilt) %just try one other filter 
        countgood = 'n'     
        StringentFilter = 'y'
        FiltersTried=FiltersTried+1
        twoImgDropParams  %reload params for next filtering
    elseif and (size(Img_Propfilt, 1)>MaxWorms*MaxWormFactor, FiltersTried == MaxFilt) 
        % if you have tried all teh filters and htere are still too many
        % worms, discard the time point
        countgood = 'n' %skip image display and mark for continue when out of loop
        break % next iteration without saving in data
    else% WORMS BELOW COUNT THRESHOLD size(Img_Propfilt, 1)<=MaxWorms*MaxWormFactor 
        countgood = 'y'  
    end   
    
        
        
end
  

if strcmpi('y', countgood)
%% collect chosen images in a matrix
  Imagesfilt={}
    for ChosenIMAGE=1:length(filtVal);
         Imagesfilt = [Imagesfilt; F(filtVal(ChosenIMAGE,1)).Image];
    end
numObj=length (Img_Propfilt (:,1))    

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
figure; imshow(imgdsp2); %axis equal;
hold on;
plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'w*')

figure ('position', scrsz); imagesc(imgdsp); %axis equal;
hold on;
plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*')

%ADD PARTICLE NUMBER
for rowN = 1:length(Img_Propfilt(:,1))
    text(Img_Propfilt(rowN,6),...
    Img_Propfilt(rowN,7),...
    num2str(Img_Propfilt(rowN,4)),'Color', 'r', 'FontSize', 12)
end

%
figure ('position', scrsz); imagesc(imgdsp); %axis equal;
hold on;
plot (Image_PropertiesAll(:,6), Image_PropertiesAll(:,7),'wO')
%plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'b*')
%ADD PARTICLE NUMBER
end


Image_PropertiesList = Image_PropertiesAll; ctX = 6; ctY = 7; TxtCol = 4; NumCol = 'r'; FntSz = 6;
if strcmpi ('y', allow_img);
    cd(CodeDir)
    PrintObjParamsAug11
end

cd([Alldata '/' TrialName 'RUNfinal'])
%% FILTERED FIGURE FOR SCORING DIAGNOSIS
imageName=[DateFldrNms2{ImN}]
if strcmpi(ProofingImgs, 'y')

%SCORED OBJECTS MULTIPLE VIEWS
figure ('position', scrsz); 
subplot(2,2,1); imshow(uint8(img1))
   hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 5)
%subplot(2,2,2); imshow(uint8(img2))
%   hold on;
%    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 5)
subplot(2,2,3); imagesc(img); title('subtracted images scaled '); colorbar
    hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 2)
subplot(2,2,4); imshow(img); title('subtracted images unscaled '); colorbar
    hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 2);...
    title (['FilterNumber_',num2str(FiltersTried)])
         
          %saveas (gcf, [DateFldrNms{W},'MultView',DateFldrNms2{ImN}(1:end-5),'.pdf'])% save as PDF
          saveas (gcf, [DateFldrNms{W},'MultView',imageName,'.pdf'])% save as PDF
   
    
%% SCORED OBJECTS MULTIPLE Parameters
%shared params
NumCol = 'K';
figure ('position', scrsz);
subplot(2,3,1); imagesc(img);
   hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r>', 'MarkerSize', 5);...
    title((['area',num2str(LowLim),'-' ,num2str(UpLim) ,'__', timeintv])); colorbar
    Image_PropertiesList = Img_Propfilt; ctX = 6; ctY = 7; TxtCol = 3
    cd(Programs)
    PrintObjParams
subplot(2,3,2); imagesc(img);
   hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r>', 'MarkerSize', 5)
    title((['MajorAxis  ',num2str(MajAxL),'-' ,num2str(MajAxU)])); colorbar
    Image_PropertiesList = Img_Propfilt; ctX = 6; ctY = 7; TxtCol = 8
    PrintObjParams  
subplot(2,3,3); imagesc(img); title('subtracted images scaled '); colorbar
    hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r.', 'MarkerSize', 5)
     title((['MinorAxis  ',num2str(MinAxL),'-' ,num2str(MinAxU)])); colorbar
    Image_PropertiesList = Img_Propfilt; ctX = 6; ctY = 7; TxtCol = 9
    PrintObjParams
subplot(2,3,4); imagesc(MasImg);
    hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'ro', 'MarkerSize', 5)
    title((['bndBx4(col14)-MasImg',num2str(BndBx4L),'-' ,num2str(BndBx4U)])); colorbar
    Image_PropertiesList = Img_Propfilt; ctX = 6; ctY = 7; TxtCol = 14
    PrintObjParams
subplot(2,3,5); imagesc(imgBW);
   hold on;
   plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'ro', 'MarkerSize', 5)
    title((['Extent-imgBW',num2str(ExtL),'-' ,num2str(ExtU)])); colorbar
    Image_PropertiesList = Img_Propfilt; ctX = 6; ctY = 7; TxtCol = 15
    PrintObjParams
subplot(2,3,6); imshow(imgdsp2); %axis equal;
    hold on;
    plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'w>')
%%    
cd([Alldata '/' TrialName 'RUNfinal'])
%saveas (gcf, [DateFldrNms{W},'ObjParams',DateFldrNms2{ImN}(1:end-5),'.pdf'])% save as PDF
saveas (gcf, [DateFldrNms{W},'ObjParams',imageName,'.pdf'])% save as PDF
%%
cd([Alldata '/' TrialName 'SubtrPics'])
figure; imagesc(img); title('subtracted images scaled '); colorbar
hold on;
plot (Img_Propfilt(:,6), Img_Propfilt(:,7),'r*', 'MarkerSize', 5)
saveas (gcf, [DateFldrNms{W},'ScaledSingle',imageName,'.pdf'])% save as PDF

end %end proofing images loop

%% CAPTURE AS STACKED TIFF
%cd([Alldata '/' TrialName 'ScoredPics'])
%Plot ontop of an image and save as a bitmapped image at the same resolution as the original image.
%This preserves the image quality and converts the vector plot objects
if strcmpi(DataCapMode, 'StackTiff')

cd([Alldata '/' TrialName 'RUNfinal'])%cd([Alldata, '/',DateFldrNms{W}]);

Nm = [strrep(DateFldrNms2{ImN}([1:end-15]), '_','-'), DateFldrNms{W}];
Ti = 'Image 1'
F = figure('Visible', proofingImgVIS);  % Create a new figure without displaying it
plotP = 'n';
plot4StacktiffAug11 (plotP, 6,7, uint8(img1), F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
I = imread('Image.tif');
    if exist([DateFldrNms{W},'stack.tif']) >0
        imwrite(I, [DateFldrNms{W},'stack.tif'], 'writemode', 'append');
    else
        imwrite(I, [DateFldrNms{W},'stack.tif']);
    end    

    
%Nm = [strrep(DateFldrNms2{ImN+1}([1:end-15]), '_','-'), DateFldrNms{W}];
%Ti = 'Image 2';
%F = figure('Visible', proofingImgVIS);  % Create a new figure without displaying it
%plotP = 'n';
%plot4StacktiffAug11 (plotP, 6,7, uint8(img2), F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
%I = imread('Image.tif');
%imwrite(I, [DateFldrNms{W},'stack.tif'], 'writemode', 'append');


%Nm = ([strrep(imageName, '_','-'),'-count',num2str(ceil(numObj/2)), '  Day ',num2str(AbsTimeDays, 3)]);
%Ti = ['Img1 + img2, SCORED', '  Day ',num2str(AbsTimeDays, 3)];
%F = figure('Visible', proofingImgVIS);  % Create a new figure without displaying it
%plotP = 'y';
%plot4StacktiffAug11 (plotP, 6,7, uint8(img1/2+img2/2), F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
%I = imread('Image.tif');
%imwrite(I, [DateFldrNms{W},'stack.tif'], 'writemode', 'append')

%Nm = ([strrep(imageName, '_','-'),'-count',num2str(ceil(numObj/2)), '  Day ',num2str(AbsTimeDays, 3)]);
%Ti = ['Img1-img2, SCORED'];
%plotP = 'y';
%F = figure('Visible', proofingImgVIS);  % Create a new figure without displaying it
%plot4StacktiffAug11 (plotP, 6,7, uint8(img), F, Nm, Ti, Img_Propfilt,'imagesc', 'Image.tif');
%I = imread('Image.tif');
%imwrite(I, [DateFldrNms{W},'stack.tif'], 'writemode', 'append')


%
Nm = ([strrep(imageName, '_','-'),'- count',num2str(ceil(numObj)),]);
Ti = ['img-', num2str(W),'/', num2str(ImN), '   TH ', num2str(Minbnd), 'to', num2str(Maxbnd), 'TH ', num2str(thresh_hold, 3)]
plotP = 'y';
F = figure('Visible', proofingImgVIS);  % Create a new figure without displaying it
plot4StacktiffAug11 (plotP, 6,7, imgdsp2, F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
hold on
%set{'Ticklabel.Color', [1,1,1]};
axes ('position', [.5, .2, .3, .3])
plot(mean(img))
text(500,500, ['\fontsize{10}', '\color{red}',num2str(min(mean(img))),'to',num2str(max(mean(img)))]) 
I = imread('Image.tif');
imwrite(I, [DateFldrNms{W},'stack.tif'], 'writemode', 'append')

end
%% EXTRACT THE CELL OF INTEREST
boundingBox=Img_Propfilt (1, 11:14)
% trim bounding box to get the proper crop
boundingBox(3)= boundingBox(3)-1
boundingBox(4)= boundingBox(4)-1
boundingBoxOffset=boundingBox
boundingBoxOffset(2)=boundingBox(2)+boundingBox(4)
%Just ofset it downwards so you avoid the intestine
%>>>boundingBoxOffset(1)=boundingBox(1)+boundingBox(3)
%% get cell intensity data
ImgCellPic=double(imcrop(img1, boundingBox));
ImgBkgdPic= double(imcrop(img1, boundingBoxOffset));
ImgCellmsk=double(Imagesfilt{1,1});
ImgCell=(ImgCellmsk.*ImgCellPic);
ImgBkgdCell=(ImgCellmsk.*ImgBkgdPic);
%% Display the croped cell...
  %>  if strcmpi (allow_img, 'y')
    figure; subplot (2,2,1); imagesc(ImgCellPic);colorbar
    subplot (2,2,2); imagesc(ImgCellmsk);colorbar
    subplot (2,2,3); imagesc(ImgCell); title('cell'); colorbar
    subplot (2,2,4); imagesc(ImgBkgdCell);title ('background'); colorbar
    
   %>I = imread('Image.tif');
    saveas(gca, ['cellhtMp', imageName,'.pdf'])
    %> input ('is this good? "y/n"')
    %> end

    
Nm = ([strrep(imageName, '_','-'),'- count',num2str(ceil(numObj)),]);
Ti = ['img-', num2str(W),'/', num2str(ImN), '   TH ', num2str(Minbnd), 'to', num2str(Maxbnd), 'TH ', num2str(thresh_hold, 3)]
plotP = 'y';
F = figure('Visible', proofingImgVIS);  % Create a new figure without displaying it
plot4StacktiffAug11 (plotP, 6,7, imgdsp2, F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
hold on
%set{'Ticklabel.Color', [1,1,1]};
axes ('position', [.5, .2, .3, .3])
plot(mean(img))
text(500,500, ['\fontsize{10}', '\color{red}',num2str(min(mean(img))),'to',num2str(max(mean(img)))]) 
saveas(gca,  ['scoredCell', imageName,'.pdf'])
 
% display cropping postions and record    
figure; imshow(img1);
hold on; 
%Imgsub=imagesc(ImgCell) %trying to plot the heatmap here
subimage(boundingBox(1), boundingBox(2), ImgCell)
subimage(boundingBoxOffset(1), boundingBoxOffset(2), ImgBkgdCell)
saveas(gca, ['SelecCell', imageName,'.pdf'])
%% extract intensity data
areaboundingBox= boundingBox(3)* boundingBox(4)
ImgCellNoZero=ImgCell(ImgCell~=0);
%figure;imshow(ImgCellNoZero);
areacell = (Img_Propfilt(1,3))
totalIntensBB = sum(sum(ImgCell))
totalIntensCell=sum(ImgCellNoZero)
MaxIntensCell = max(ImgCellNoZero)
MinIntensCell = min(ImgCellNoZero)
rangecell=range(ImgCellNoZero)
stdevcell=std2(ImgCellNoZero)
AvgIntensCell = totalIntensCell/areacell 
totalIntensBkgd = sum(sum(ImgBkgdCell))
AvgIntensBkgd= totalIntensBkgd/areacell 

totalIntensAbvBkgd = totalIntensCell-totalIntensBkgd;
AvgIntensAbvBkgd = AvgIntensCell-AvgIntensBkgd
%MaxIntensAbvBkg = max(max(totalIntensAbvBkgd))
%%
%>I = imread('Image.tif');
%>imwrite(I, [DateFldrNms{W},'stack.tif'], 'writemode', 'append')
%% GetSUBIMAGE OF CELL
%% save compile data
cd([Alldata '/' TrialName 'RUNfinal'])
save ([imageName, '.mat'], 'img1','Image_PropertiesAll',...
    'Image_PropertiesList', 'Images', 'Imagesfilt', 'Img_Propfilt') %, 'CumulImg' 
%figure; imagesc(img)
end

close all
