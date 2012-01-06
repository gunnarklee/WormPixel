function [Totalcelldata cellnums] = OneWormAnalysis(Imagesfilt, Img_Propfilt, img1)
%UNTITLED2 Summary of this function goes here
%Single worm is input, get image characteristics
pad=20


if iscell(Imagesfilt)
WmImg=Imagesfilt{1,1};
else 
WmImg=Imagesfilt;
end

%WmImg=~WmImg
%figure ;imshow(WmImg)
size(WmImg)%pad image
WmImgPad = zeros(size(WmImg)+pad)
Rowpad=(pad/2+1:size(WmImgPad, 1)-pad/2)
ColPad=(pad/2+1:size(WmImgPad, 2)-pad/2)
WmImgPad(Rowpad, ColPad)=WmImg
%figure ;imshow(WmImgPad)
%%   get spine
close all
SE=strel('disk', 3);
minWm=imerode(WmImgPad, SE);

figure ;imshow(minWm, 'InitialMagnification', 400); title ('minworm')
figure ;imshow(WmImgPad,'InitialMagnification', 400); title ('WmImgPad')

%%
[x,y]=find(minWm)
xx=1:size(WmImgPad, 2)
size(WmImgPad)


skele=bwmorph(WmImgPad, 'skel', Inf);
%[x,y]=ind2sub(size(skele), find(skele))
figure; imshow(imoverlay (mat2gray(WmImgPad), skele,  [255, 0, 0]));
skele2=bwmorph(skele, 'spur');
figure; imshow(imoverlay (mat2gray(WmImgPad), skele2,  [0,255, 0]));

skele3=bwmorph(skele2, 'spur');
figure; imshow(imoverlay (mat2gray(WmImgPad), skele3,  [0, 0, 255]));

skeleSH=bwmorph(skele3, 'shrink');
figure; imshow(imoverlay (mat2gray(WmImgPad), skeleSH,  [0, 0, 255]), 'InitialMagnification', 400)
[x,y]=ind2sub(size(skeleEND), find(skeleSH))

%% build list of endpoints

skeleEND=bwmorph(skeleSH, 'endpoints');
figure; imshow(imoverlay (mat2gray(WmImgPad), skeleEND,  [0, 0, 255]), 'InitialMagnification', 400)
[x,y]=ind2sub(size(skeleEND), find(skeleEND))
seed= [x(1),y(1)]




%figure; imshow(WmImgPad)
%[Totalcelldata cellnums]=cloadfinal2(WmImgPad)
%figure; imshow(ims)
%hold on
%[Totalcelldata cellnums]=cload(WmImgPad)
%
%%sicic spline finding and segmenting
cplot(Totalcelldata, 'R', 'mid', 1)
cplot(Totalcelldata, 'b', 'border', 1)
%cplot(celldata, str, feat,cellnums,  width)


%% extracting original image
%get bounding box
Xmin=Img_Propfilt(11)
Xmax=Img_Propfilt(11)+Img_Propfilt(13)
Ymin=Img_Propfilt(12)
Ymax=Img_Propfilt(12)+Img_Propfilt(14)
% get worm image
%figure; imshow(img1(Ymin:Ymax, Xmin:Xmax))
ims=img1(Ymin:Ymax, Xmin:Xmax)


end

