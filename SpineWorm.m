function [SpineData, poshead] = SpineWorm (WmImgPad, img1, allow_img, stoppoint, poshead, numpts)
%Single worm is input, get SPINE characteristics

%% GET SPINE
xx=1:size(WmImgPad, 2)
size(WmImgPad)

skele=bwmorph(WmImgPad, 'skel', Inf);
skele2=bwmorph(skele, 'spur');
[DiffPts1]= MtxDiff (skele, skele2); %record endpoints removed
skele3=bwmorph(skele2, 'spur');
[DiffPts2]= MtxDiff (skele2, skele3); %record endpoints removed
skeleSH=bwmorph(skele3, 'shrink');
[DiffPts3]= MtxDiff (skele3, skeleSH); %record endpoints removed
DiffPts=[DiffPts1; DiffPts2; DiffPts3]

%reconsitute endpoints
%%>> IN PROGRESS <<<<

row=1:length(DiffPts) % LIST OF ENDPOINTS REMOVED
row=sort(row, 'descend')
%>>for PointTst=row
%>>DiffPts(PointTst)   
%>> skeleCurr 
%>> END=bwmorph(skeleCurr, 'endpoints')
%>>end
%%>> IN PROGRESS <<<<



if (strcmpi (allow_img, 'y'));  
  %  figure ;imshow(minWm, 'InitialMagnification', 400); title ('minworm')
    figure ;imshow(WmImgPad,'InitialMagnification', 400); title ('WmImgPad')
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
    
  %>  if isempty(poshead)
  %>      WmImg=imoverlay (mat2gray(WmImgPad), WorkSpine,  [0, 0, 255]);title ('skeleSH-shrunk');
  %>      [poshead] = GetPoint(WmImg, [ceil(size(WmImg (:,:,1), 2)*.85), ceil(size(WmImg (:,:,1), 1)*.85)])
  %>  end
    
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
    
    %>>>>>>>SpineData.AngleLs=GetAngles(Pointlist)
    
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
    SpineData.Lengthls = Lengthls
    %SpineData.PtpointdistLs =PtpointdistLs
    SpineData.SpineList = SpineList
    SpineData.Pointlist = Pointlist
    SpineData.endpoints = size(x, 1)
    SpineData.poshead=poshead
    SpineData.padimg=WmImgPad
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
