function SpineWorm(Imagesfilt, Img_Propfilt, img1)
%UNTITLED2 Summary of this function goes here
%Single worm is input, get image characteristics
pad=20
%CmapName='copper'
numpts=13 %number of spine points

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
figure; imshow(imoverlay (mat2gray(WmImgPad), skele,  [255, 0, 0]), 'InitialMagnification', 400); title ('skele-original');
skele2=bwmorph(skele, 'spur');
figure; imshow(imoverlay (mat2gray(WmImgPad), skele2,  [0,255, 0]), 'InitialMagnification', 400); title ('skele2-despurred1x');

skele3=bwmorph(skele2, 'spur');
figure; imshow(imoverlay (mat2gray(WmImgPad), skele3,  [0, 0, 255]), 'InitialMagnification', 400); title ('skele3-despurred2x');

skeleSH=bwmorph(skele3, 'shrink');
figure; imshow(imoverlay (mat2gray(WmImgPad), skeleSH,  [0, 0, 255]), 'InitialMagnification', 400);title ('skeleSH-shrunk');
%Losing worm ends
%May need to minimize srinling or add back ends to get full worm
%% Order the spine points
SpineList=[]

%select spine to use
WorkSpine=skeleSH;
[wkSpnX,wkspnY]=ind2sub(size(WorkSpine), find(WorkSpine));
cmap=colormap(jet (size(wkSpnX,1))) %copper
WmImgPadcolor=mat2gray(WmImgPad)
Lengthls=[]

for Pt=1:size(wkSpnX, 1)-1
    %close all
    skeleEND=bwmorph(WorkSpine, 'endpoints');
    [x,y]=ind2sub(size(skeleEND), find(skeleEND));
    CurrPt= [x(1),y(1)];
    anchor =[x(2),y(2)]; %find anchor and new endpoint
    
    
    if size (SpineList, 1) > 0
        lastPt=SpineList(end, :) %get last point appended
        %if x and y values for both points are different get
        %you have a slanted line, find the hypotenuse
        if (CurrPt(1) == lastPt(1)) + (CurrPt(2) == lastPt(2)) == 0
            ThridPt=[CurrPt(1),lastPt(2)]
            
            %This works because points are in a line with the third point
            distA=abs(CurrPt(1)-ThridPt(1)+CurrPt(2)-ThridPt(2))
            distB=abs(lastPt(1)-ThridPt(1)+lastPt(2)-ThridPt(2))
            dist=hypot(distA, distB)
        else % the points are direclty next to each other
            dist=1
        end
        %if
        %% if both the x and y values are differnt get hypotenuse
        %distA=((lastPt(1,1)-CurrPt(1,1))^2+(lastPt(1,2)-CurrPt(1,2))^2)
        %distb=((lastPt(1,1)-CurrPt(1,1))^2+(lastPt(1,2)-CurrPt(1,2))^2)
        
        %dist=sqrt((lastPt(1,1)-CurrPt(1,1))^2+(lastPt(1,2)-CurrPt(1,2))^2)
        % if only one coordinate value is diff, get the linear dist
        
        Lengthls=[Lengthls; dist+Lengthls(end)] %get cumulative length
    else
        Lengthls=0
    end
    
    
    %append the current point
    SpineList=[SpineList;CurrPt];
    
    %addpoint to imgage with new color
    WmImgPadcolor=(imoverlay (mat2gray(WmImgPadcolor), skeleEND,  cmap(Pt,:)));
    
    %erase the current point
    WorkSpine(CurrPt(1),CurrPt(2))=0;
end

%Finsh with last point
SpineList=[SpineList;anchor];
%figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);

%% get distances and N points along spine

SegmtLn=Lengthls(end)/numpts
pointloc=[SpineList(1,:)] %start at the beginning
PtDistLs=[]

for SpPt=1:numpts
    %SpPt=1
    PtDist=SpPt*SegmtLn
    %capture point coordinates at the point before the segment exceedes the
    %desired size
    pointloc=[pointloc; SpineList(min(find(Lengthls >PtDist))-1, :)]
    PtDistLs=[PtDistLs, PtDist]
end
%capture the final point.
pointloc=[pointloc;SpineList(end,:)]
%%

%WmImgPadcolor=(imoverlay (WmImgPadcolor, pointloc,  [0,0,255]));
figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);
hold on
plot(pointloc(:,2), pointloc(:,1), 'b*', 'MarkerSize', 10)


%% use the law of cosines to get angles between points
AngleLs=[]
for Pair=1:size(PtDistLs, 2)-2

triPt=[pointloc(Pair,2), pointloc(Pair,1);...
    pointloc(Pair+1,2), pointloc(Pair+1,1);...
    pointloc(Pair+2,2), pointloc(Pair+2,1)]

%get lengths of lines between points use the distance formula
Aln = sqrt((triPt(2,1) - triPt(1,1))^2 + (triPt(2,2) - triPt(2,1))^2)
Bln = sqrt((triPt(3,1) - triPt(2,1))^2 + (triPt(3,2) - triPt(2,1))^2)
Cln = sqrt((triPt(3,1) - triPt(1,1))^2 + (triPt(3,2) - triPt(1,1))^2)

angle=acosd((Aln^2+Bln^2-Cln^2)/(2*Aln*Bln))
%angle=57.2957795*angle%in degres
AngleLs=[AngleLs;angle]
end
%% get curvature matrix and curvature column



%% extracting original image
%get bounding box
%Xmin=Img_Propfilt(11)
%%Xmax=Img_Propfilt(11)+Img_Propfilt(13)
%Ymin=Img_Propfilt(12)
%Ymax=Img_Propfilt(12)+Img_Propfilt(14)
% get worm image
%figure; imshow(img1(Ymin:Ymax, Xmin:Xmax))
%ims=img1(Ymin:Ymax, Xmin:Xmax)


collect important outputs 


end

