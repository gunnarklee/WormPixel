function [AngleLs]=GetAngles(Pointlist, allow_img)

%% Pierce-Shimomoura's "curvature column"
    %% use the law of cosines to get angles between points
    load 
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
        AngleLs=[AngleLs;angle]; %,angleB, AnglDir]
        
        %% DIAGNOSITICS
        if strcmpi(allow_img, 'y')
            
            %get lengths of lines between points use the pointdistance formula
            ABlength = sqrt((PointA(1) - PointB(1))^2 + (PointA(2) - PointB(2))^2);
            BClength = sqrt((PointB(1) - PointC(1))^2 + (PointB(2) - PointC(2))^2);
            CAlength = sqrt((PointC(1) - PointA(1))^2 + (PointC(2) - PointA(2))^2);
            
            
            %Change new coordiantes to cartesain FOR PLOTTING and checking
            [CnewPt(1),CnewPt(2)]=pol2cart(CnewPolar(1),CnewPolar(2));
            [AnewPt(1),AnewPt(2)]=pol2cart(AnewPolar(1),AnewPolar(2));
            
            figure; imshow(img1)
            hold on
            plot(PointAori(1),PointAori(2),'b+',  PointBori(1), PointBori(2), 'b*',...
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
            
            
            
%% DIAGNOSTIC- check for consistencies in old versus new triangles -
            if strcmpi (stoppoint, 'y')
            %TRIG VERIFICATION
            angleB=acosd((ABlength^2+BClength^2-CAlength^2)/(2*ABlength*BClength));
            ABNewlength = sqrt((AnewPt(1) - PointBori(1))^2 + (AnewPt(2) - PointBori(2))^2);
            BCNewlength = sqrt((PointBori(1) - CnewPt(1))^2 + (PointBori(2) - CnewPt(2))^2);
            CANewlength = sqrt((CnewPt(1) - AnewPt(1))^2 + (CnewPt(2) - AnewPt(2))^2);
            angleBnew=acosd((ABNewlength^2+BCNewlength^2-CANewlength^2)/(2*ABNewlength*BCNewlength));
            oldvsnew=[ABlength,ABNewlength;BClength,BCNewlength;CAlength,CANewlength;angleB,angleBnew];
            %x = xCenter + radius * cos(angle)
            %y = yCenter + radius * sin(angle)
            
                stopPt= input ('next image?', 's')
            end
        end %end DIAGNOSTIC SECTION
    end
    
end    