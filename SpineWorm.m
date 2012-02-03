function [SpineData, poshead] = SpineWorm(Imagesfilt, Img_Propfilt, img1, ErrorDir, allow_img, stoppoint, poshead)
%UNTITLED2 Summary of this function goes here
%Single worm is input, get image characteristics
pad=20
numpts=12 %number of spine points

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
skele2=bwmorph(skele, 'spur');
skele3=bwmorph(skele2, 'spur');
skeleSH=bwmorph(skele3, 'shrink');

if (strcmpi (allow_img, 'y'));
    figure; imshow(imoverlay (mat2gray(WmImgPad), skele,  [255, 0, 0]), 'InitialMagnification', 400); title ('skele-original');
    figure; imshow(imoverlay (mat2gray(WmImgPad), skele2,  [0,255, 0]), 'InitialMagnification', 400); title ('skele2-despurred1x');
    figure; imshow(imoverlay (mat2gray(WmImgPad), skele3,  [0, 0, 255]), 'InitialMagnification', 400); title ('skele3-despurred2x');
    figure; imshow(imoverlay (mat2gray(WmImgPad), skeleSH,  [0, 0, 255]), 'InitialMagnification', 400);title ('skeleSH-shrunk');
end
%Losing worm ends
%May need to minimize srinking or add back ends to get full worm


%% Endpoint ERRORS - fix or catch them
skeleEND=bwmorph(skeleSH, 'endpoints');
[x,y]=ind2sub(size(skeleEND), find(skeleEND));

if (size (x, 1) == 2) == 0 %if it does not equal 2 as expected
    if (size (x, 1) > 2) == 1 % The spine has spurs
        %spurfix  *********
    end
    
    if (size (x, 1) < 2) == 1 %less than two endoints worm is a circle
        %find thinnest point on worm body and delete the connecting pixel
        %******** threeconnected pint deletion ?
    end
    
    %if the fixes did not create two endpoints save particle for later
    if (size (x, 1) == 2) == 0; SpineData.spinegood ='n'; end
else
    SpineData.spinegood ='y'
end


%% Order the spine points
if strcmpi (SpineData.spinegood, 'n')  == 0 % if the spine is good, proceede
    SpineList=[]
    Distlist=[]
    %select spine to use
    WorkSpine=skeleSH;
    [wkSpnX,wkspnY]=ind2sub(size(WorkSpine), find(WorkSpine));
    cmap=colormap(jet (size(wkSpnX,1))) %copper
    WmImgPadcolor=mat2gray(WmImgPad)
    Lengthls=[]
    
    
    %Get the head point and assign anterior, posterior
    %if point exists from last image use that point
    %for Pt=1:size(SpineData.SpineList, 1)-1
    %figure; imoverlay(Imagesfilt{1,1}, [SpineData.SpineList(Pt,1),SpineData.SpineList(Pt,2)], cmap(Pt,:));
    %hold on
    %end
    
    if isempty(poshead)
        WmImg=imoverlay (mat2gray(WmImgPad), WorkSpine,  [0, 0, 255]);title ('skeleSH-shrunk');
        [poshead] = GetPoint(WmImg, [ceil(size(WmImg (:,:,1), 2)*.85), ceil(size(WmImg (:,:,1), 1)*.85)])
    end
    
    %set anchor once use ~(not)anchor as curr point
    skeleEND=bwmorph(WorkSpine, 'endpoints');
    [x,y]=ind2sub(size(skeleEND), find(skeleEND));
    

    %Distance between head position and top of matrix
    TopHeaddist=sqrt((poshead(2) - x(1))^2 + (poshead(1) - y(1))^2)
    %Distance between head position and bottom of matrix
    BotHeaddist=sqrt((poshead(2) - x(2))^2 + (poshead(1) - y(2))^2)
    %reorder top of matrix to correspond with known head postion
    
    
    %choose the point furthest from the head to be the anchor
    if BotHeaddist<TopHeaddist % then the head is at the bottom of the matrix
        anchor =[x(1),y(1)]; % flip matix my making the first endpoint into the last
    else
        anchor =[x(2),y(2)]; % the matirx is in the correct orientation, anchor the bottom
    end
    
%% BUILD SPINE - use ~(not)anchor point as current point    
    for Pt=1:size(wkSpnX, 1)-1
        skeleEND=bwmorph(WorkSpine, 'endpoints');
        [x,y]=ind2sub(size(skeleEND), find(skeleEND));
        
        % Identify the position of the anchor row
        anchorROW=find (x == anchor(1))
        if length(anchorROW) > 1 % if both endpoints have the same X
            anchorROW=find (y == anchor(2))
        end
        
    %for each iteration of spine decay check that there are still 2 endpoints
    %deals with case of spine that crosses itself to make a circle
    if size (x,1) <2  % bad spine case
    SpineData.spinegood = 'n'
    SpineData.endpoints = size(x, 1);
    break
    end
    
        
        
        % Take the non-anchor point as the current point
        % if it can't locate the anchor at any iteration the spine is bad
        
        if anchorROW == 1
            CurrPt= [x(2),y(2)];
        else
            CurrPt= [x(1),y(1)];
        end
        
        if size (SpineList, 1) > 0
            lastPt=SpineList(end, :) %get last point appended
            %if x and y values for both points are different get
            %you have a slanted line, find the hypotenuse
            %             if (CurrPt(1) == lastPt(1)) + (CurrPt(2) == lastPt(2)) == 0
            %                 ThridPt=[CurrPt(1),lastPt(2)]
            %
            %                 %This works because points are in a line with the third point
            %                 pointdistA=abs(CurrPt(1)-ThridPt(1)+CurrPt(2)-ThridPt(2))
            %                 pointdistB=abs(lastPt(1)-ThridPt(1)+lastPt(2)-ThridPt(2))
            %                 pointdist=hypot(pointdistA, pointdistB)
            %             else % the points are direclty next to each other
            %                 pointdist=1
            %             end
            
            %get the distance between the two points
            pointdist=sqrt((lastPt(1) - CurrPt(1))^2 + (lastPt(2) - CurrPt(2))^2)
            Lengthls=[Lengthls; pointdist+Lengthls(end)] %get cumulative length
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
    pointdist = sqrt((anchor(1) - CurrPt(1))^2 + (anchor(2) - CurrPt(2))^2)
    Lengthls=[Lengthls; pointdist+Lengthls(end)]
    SpineList=[SpineList;anchor];
    
    %[SpineList, Lengthls]
    
    %figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);
    
    %% get pointdistances and N points along spine
    
    SegmtLn=Lengthls(end)/(numpts-1)
    Pointlist=[SpineList(1,:)] %start at the beginning
    %Pointlist=[pointloc]
    Distlist=[Distlist;0, Lengthls(1, :) ]
    for SpPt=1:numpts-2 %-2,  first & last points added outside loop
        %SpPt=1
        Ptpointdist=SpPt*SegmtLn
        %get point coord at the point before the segment exceedes the
        %desired size
        RowHit=min(find(Lengthls >Ptpointdist))-1
        Pointlist=[Pointlist; SpineList(RowHit, :)]
        Distlist=[Distlist;Ptpointdist, Lengthls(RowHit, :) ]% list target and realized distances
    end
    %capture the final point.
    Pointlist=[Pointlist;SpineList(end,:)]
    Distlist=[Distlist;(numpts-1)*SegmtLn, Lengthls(end, :) ]% list target and realized distances
    
    %%
    
    %WmImgPadcolor=(imoverlay (WmImgPadcolor, pointloc,  [0,0,255]));
    if strcmpi (allow_img, 'y')
        figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);
        hold on
        plot(Pointlist(:,2), Pointlist(:,1), 'b*', 'MarkerSize', 10)
    end
    
    
    %% Pierce-Shimomoura's "curvature column"
    %% use the law of cosines to get angles between points
    
    AngleLs=[]
    for Pair=1:size(Pointlist, 1)-2
        
        %get three point locations A near end, C nearest anchor
        PointA=[Pointlist(Pair,1), Pointlist(Pair,2)];
        PointB=[Pointlist(Pair+1,1), Pointlist(Pair+1,2)];
        PointC=[Pointlist(Pair+2,1), Pointlist(Pair+2,2)];
        
        
        %% for each three points get angle and side relative to bearing
        %linear trnsform angle so point B is at origin
        PointBori=[0,0];
        PointAori=[PointA(1)-PointB(1),PointA(2)-PointB(2)];
        PointCori=[PointC(1)-PointB(1),PointC(2)-PointB(2)];
        
        %transform orig. points to the cartesion axis
        [PtATheta,PtARoh]=cart2pol(PointAori(1),PointAori(2));
        %       [PtATheta] =MkClockWise(PtATheta)
        [PtCTheta,PtCRoh]=cart2pol(PointCori(1),PointCori(2));
        %       [PtCTheta] =MkClockWise(PtCTheta)
        
        %find angle required to rotate point C to the left of B (to 180deg)
        [MoveCto180Theta,MoveCto180_Roh]=cart2pol(-PtCRoh, 0); % position th the left of
        %      [PosLeftBTheta] =MkClockWise(PosLeftBTheta)
        
        CRotation=MoveCto180Theta-PtCTheta; %rotn from C_orig to get to the 180 mark
        CnewPolar=[MoveCto180Theta,MoveCto180_Roh]; %rotate C to the 180
        AnewPolar=[CRotation+PtATheta,PtARoh]; %rotate A by the same amount
        
        %determine polarity of inflection of angle ABC
        %      if AnewPt(1) > 0, AnglDir = -1, else AnglDir = 1, end
        % get central angle and check
        
        angle=AnewPolar(1)*(180/pi)%%-360 % get angle and assign polarity
        if abs(angle) > 180 % constrain to under 180 degrees
            angle=angle-360
        end    
        AngleLs=[AngleLs;angle] %,angleB, AnglDir]
        
        %% DIAGNOSITICS
        if strcmpi(allow_img, 'y')
            
            %get lengths of lines between points use the pointdistance formula
            ABlength = sqrt((PointA(1) - PointB(1))^2 + (PointA(2) - PointB(2))^2);
            BClength = sqrt((PointB(1) - PointC(1))^2 + (PointB(2) - PointC(2))^2);
            CAlength = sqrt((PointC(1) - PointA(1))^2 + (PointC(2) - PointA(2))^2);
            
            
            %Change new coordiantes to cartesain FOR PLOTTING and checking
            [CnewPt(1),CnewPt(2)]=pol2cart(CnewPolar(1),CnewPolar(2));
            [AnewPt(1),AnewPt(2)]=pol2cart(AnewPolar(1),AnewPolar(2));
            
            
            figure; plot(PointAori(1),PointAori(2),'b+',  PointBori(1), PointBori(2), 'b*',...
                PointCori(1), PointCori(2), 'bo')
            hold on
            plot(PointA(1),PointA(2),'r+',  PointB(1), PointB(2), 'r*',...
                PointC(1), PointC(2), 'ro')
            hold on
            plot(AnewPt(1),AnewPt(2),'g+',  PointBori(1), PointBori(2), 'g*',...
                CnewPt(1),CnewPt(2), 'go')
            
            
            text(PointA(1),PointA(2), 'A');
            text(PointB(1),PointB(2), 'B');
            text(PointC(1),PointC(2), 'C');
            %text(PointCNew(1),PointCNew(2), 'CNew');
            %text(PointANew(1),PointANew(2), 'ANew');
            
            text(AnewPt(1),AnewPt(2), 'A');
            text(0,0, 'B');
            text(CnewPt(1),CnewPt(2), 'C');
            
            text(PointAori(1),PointAori(2), 'A');
            text(PointBori(1),PointBori(2), 'B');
            text(PointCori(1),PointCori(2), 'C');
            text(PointBori(1)-10,PointBori(2)+15, ['angle', num2str(angle)]);
            
            xlim([-20,100]);
            ylim([-20,100]);
            
            
            
            %DIAGNOSTIC- check for consistencies in old versus new triangles -
            %TRIG VERIFICATION
            angleB=acosd((ABlength^2+BClength^2-CAlength^2)/(2*ABlength*BClength));
            ABNewlength = sqrt((AnewPt(1) - PointBori(1))^2 + (AnewPt(2) - PointBori(2))^2);
            BCNewlength = sqrt((PointBori(1) - CnewPt(1))^2 + (PointBori(2) - CnewPt(2))^2);
            CANewlength = sqrt((CnewPt(1) - AnewPt(1))^2 + (CnewPt(2) - AnewPt(2))^2);
            angleBnew=acosd((ABNewlength^2+BCNewlength^2-CANewlength^2)/(2*ABNewlength*BCNewlength));
            oldvsnew=[ABlength,ABNewlength;BClength,BCNewlength;CAlength,CANewlength;angleB,angleBnew];
            %x = xCenter + radius * cos(angle)
            %y = yCenter + radius * sin(angle)
            if strcmpi (stoppoint, 'y')
                stopPt= input ('next image?', 's')
            end
        end %end DIAGNOSTIC SECTION
    end
    
    %get bounding box
    %Xmin=Img_Propfilt(11)
    %%Xmax=Img_Propfilt(11)+Img_Propfilt(13)
    %Ymin=Img_Propfilt(12)
    %Ymax=Img_Propfilt(12)+Img_Propfilt(14)
    % get worm image
    %figure; imshow(img1(Ymin:Ymax, Xmin:Xmax))
    %ims=img1(Ymin:Ymax, Xmin:Xmax)
    
    
    %%collect important outputs
    SpineData.spinegood = 'y'
    SpineData.AngleLs = AngleLs;
    SpineData.Lengthls = Lengthls
    %SpineData.PtpointdistLs =PtpointdistLs
    SpineData.SpineList = SpineList
    SpineData.Pointlist = Pointlist
    SpineData.endpoints = size(x, 1)
else  % end good spine process
    
    % bad spine case
    SpineData.spinegood = 'n'
    SpineData.endpoints = size(x, 1);
end
end


function [CWangle] =MkClockWise(Agl)
%take in a polar coordinate angle, check rotation direction, make clockwise
if Agl < 0,
    CWangle=(2*(pi) + Agl)
else
    CWangle=Agl
end

end
