function [SpineData, poshead] = SpineWorm (WmImgPad,allow_img, poshead, numpts)
%Single worm is input, get SPINE characteristics

%% GET SPINE
xx=1:size(WmImgPad, 2);

%%
WmImgPad=imfill(WmImgPad, 'holes');
skele=bwmorph(WmImgPad, 'skel', Inf);
%figure; imshow(imoverlay (mat2gray(WmImgPad), skele,  [255, 0, 0]), 'InitialMagnification', 400); title ('skele-original');
%%
skeleSH=bwmorph(skele, 'thin', Inf);
%%
skeleSH=bwmorph(skele, 'shrink');
endpoints = bwmorph(skeleSH, 'endpoints');
while sum(endpoints(:)) > 2
    skeleSH=bwmorph(skele, 'shrink');
    sum(endpoints(:)); %remove  will report
    skele=bwmorph(skeleSH, 'spur');
    endpoints = bwmorph(skeleSH, 'endpoints');
end


%%
%>>skele3=bwmorph(skele2, 'spur');
%>>[DiffPts2]= MtxDiff (skele2, skele3); %record endpoints removed

%>>[DiffPts3]= MtxDiff (skele3, skeleSH); %record endpoints removed
%>>DiffPts=[DiffPts1; DiffPts2; DiffPts3]

%reconsitute endpoints
%%>> IN PROGRESS <<<<

%>>row=1:length(DiffPts) % LIST OF ENDPOINTS REMOVED
%>>row=sort(row, 'descend')
%>>for PointTst=row
%>>DiffPts(PointTst)
%>> skeleCurr
%>> END=bwmorph(skeleCurr, 'endpoints')
%>>end
%%>> IN PROGRESS <<<<

%%

if (strcmpi (allow_img, 'y'));
    figure ;imshow(WmImgPad,'InitialMagnification', 400); title ('WmImgPad')
    figure; imshow(imoverlay (mat2gray(WmImgPad), skele,  [255, 0, 0]), 'InitialMagnification', 400); title ('skele-original');
    figure; imshow(imoverlay (mat2gray(WmImgPad), skeleSH,  [0, 0, 255]), 'InitialMagnification', 400);title ('skeleSH-shrunk');
end
%Losing worm ends

%% Endpoint ERRORS - fix or catch them
skeleEND=bwmorph(skeleSH, 'endpoints');
[x,y]=ind2sub(size(skeleEND), find(skeleEND));

if (size (x, 1) == 2) == 0 %if it does not equal 2 as expected
    if (size (x, 1) > 2) == 1 % The spine has spurs
        %spurfix  *********
    end
    
    if (size (x, 1) < 2) == 1 %less than two endoints worm is a circle
        %find thinnest point on worm body and delete the connecting pixel
        %******** threeconnected point deletion ?
    end
    
    %if the fixes did not create two endpoints save particle for later
    %abnd mark as bad spine
    if (size (x, 1) == 2) == 0; SpineData.spinegood ='n'; end
    SpineData.FailPt= 'endpoints';
else
    SpineData.spinegood ='y'; % if there are no errors the spine is good
end


%% Order the spine points
if strcmpi (SpineData.spinegood, 'n') == 0 % if the spine is good, proceede
    SpineList=[];
    %Distlist=[];
    %select spine to use
    WorkSpine=skeleSH;
    [wkSpnX,wkspnY]=ind2sub(size(WorkSpine), find(WorkSpine));
    
    if strcmpi (allow_img, 'y')
        cmap=colormap(jet (size(wkSpnX,1))); %copper
        WmImgPadcolor=mat2gray(WmImgPad);
    end
    
    Lengthls=[];
    
    
    %% Get the head point and assign anterior, posterior
    %if point exists from last image use that point
    
    %set anchor once use ~(not)anchor as curr point
    skeleEND=bwmorph(WorkSpine, 'endpoints');
    [x,y]=ind2sub(size(skeleEND), find(skeleEND));
    
    
    %Distance between head position and top of matrix
    TopHeaddist=sqrt((poshead(2) - x(1))^2 + (poshead(1) - y(1))^2);
    %Distance between head position and bottom of matrix
    BotHeaddist=sqrt((poshead(2) - x(2))^2 + (poshead(1) - y(2))^2);
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
        anchorROW=find (x == anchor(1));
        if length(anchorROW) > 1 % if both endpoints have the same X
            anchorROW=find (y == anchor(2));
        end
        
        %for each iteration of spine decay check that there are still 2 endpoints
        %deals with case of spine that crosses itself to make a circle
        if size (x,1) <2  % bad spine case
            SpineData.spinegood = 'n';
            SpineData.endpoints = size(x, 1);
            return
        end
        
        % Take the non-anchor point as the current point
        % if it can't locate the anchor at any iteration the spine is bad
        
        if anchorROW == 1
            CurrPt= [x(2),y(2)];
        else
            CurrPt= [x(1),y(1)];
        end
        
        if size (SpineList, 1) > 0
            lastPt=SpineList(end, :); %get last point appended
            
            %get the distance between the two points
            pointdist=sqrt((lastPt(1) - CurrPt(1))^2 + (lastPt(2) - CurrPt(2))^2);
            Lengthls=[Lengthls; pointdist+Lengthls(end)]; %get cumulative length
        else
            Lengthls=0;
        end
        
        %append the current point
        SpineList=[SpineList;CurrPt];
        
        if strcmpi (allow_img, 'y')
            %addpoint to imgage with new color
            WmImgPadcolor=(imoverlay (mat2gray(WmImgPadcolor), skeleEND,  cmap(Pt,:)));
            %erase the current point
        end
        
        WorkSpine(CurrPt(1),CurrPt(2))=0;
    end
    
    %Finsh with last point
    pointdist = sqrt((anchor(1) - CurrPt(1))^2 + (anchor(2) - CurrPt(2))^2);
    Lengthls=[Lengthls; pointdist+Lengthls(end)];
    SpineList=[SpineList;anchor];
    
    
    %% get pointdistances and N points along spine
    avgwin=3; %rounds decimals down so flank for 2=1  & flank for 3=1
    flank=floor(avgwin/2);
    
    %ERROR CHECK - short spines will not work dump to "spine error"
    if numpts*avgwin > length(SpineList)
        SpineData.spinegood = 'n';
        SpineData.endpoints = size(x, 1);
        return
    end
    
    %preallocate matricies
    Distlist=zeros(numpts,2);
    Pointlist=zeros(numpts,2);
    
    SegmtLn=Lengthls(end)/(numpts-1);
    Pointlist(1,:)=[SpineList(1,:)]; %start at the beginning
    %Pointlist=[pointloc]
    Distlist(1,:)=[0, Lengthls(1, :)];
    %Determine window size for point averaging
    %win=(ceil((size(SpineList, 1)/(size(Pointlist,1)-1))*2-1)/2)
    
    for SpPt=1:numpts-2 %-2,  first & last points added outside loop
        Ptpointdist=SpPt*SegmtLn;
        %get point coord at the point before the segment exceedes the
        %desired size
        RowHit=find(Lengthls >Ptpointdist, 1 )-1;
        
        %reduce wiggle with sliding window for point
        try %this step fails sometime
            PointLoc=median(SpineList((RowHit-flank:RowHit+flank), :));
        catch
            SpineData.spinegood = 'n';
            SpineData.endpoints = size(x, 1);
            SpineData.FailPt= 'pointloc';
            return%save('SpineFailLn179')
        end
        
        try
            Pointlist(SpPt+1,:)=PointLoc; %concatente
        catch
            SpineData.spinegood = 'n';
            SpineData.endpoints = size(x, 1);
            SpineData.FailPt= 'pointlist';
            return%save('SpineFailLn179')
        end
        
        Distlist(SpPt+1,:)=[Ptpointdist, Lengthls(RowHit, :) ];% list target and realized distances
    end
    %capture the final point.
    Pointlist(end,:)=SpineList(end,:);
    Distlist(end, :)=[(numpts-1)*SegmtLn, Lengthls(end, :) ];% list target and realized distances
    
    %% IMAGE verification
    %WmImgPadcolor=(imoverlay (WmImgPadcolor, pointloc,  [0,0,255]));
    if strcmpi (allow_img, 'y')
        figure; imshow(WmImgPadcolor, 'InitialMagnification', 800);
        hold on
        plot(Pointlist(:,2), Pointlist(:,1), 'b+', 'MarkerSize', 20);
        plot (Pointlist(:,2),Pointlist(:,1), 'b-');
        hold off
    end
    
    
    %%collect important outputs
    SpineData.spinegood = 'y';
    SpineData.Lengthls = Lengthls;
    SpineData.SpineList = SpineList;
    SpineData.Pointlist = Pointlist;
    SpineData.endpoints = size(x, 1);
    SpineData.poshead=poshead;
    SpineData.padimg=WmImgPad;
else   %bad spine case
    %end good spine process and mark as bad spine
    SpineData.spinegood = 'n';
    SpineData.endpoints = size(x, 1);
end
end


function [CWangle] =MkClockWise(Agl)
%take in a polar coordinate angle, check rotation direction, make clockwise
if Agl < 0,
    CWangle=(2*(pi) + Agl);
else
    CWangle=Agl;
end

end
