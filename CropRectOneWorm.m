%CropCircleAdptAuto3 Gunnar Kleemann (6/19/08)
%version 3 adds optional image supression
%Circular masking to select 4 regions


%check with markup if defaults are OK
lng1=length(imgtmp(1,:));
lng2=length(imgtmp(:,1));


%% 
if exist('xctr') == 0 % if the plate location is not known ask the user
%figure; imagesc(imgtmp)
GoodCrop = 'n'; % for the first iteration you will allow user to recrop until they are satissfied
% this is reqired because the masking circle is biggger than expected
while strcmpi(GoodCrop ,'n')

disp('indicate the included region');
disp('DOUBLE CLICK WHEN DONE!');
%% 
        clear ('posctr'); clear ('posedge');

        while (exist('posctr'))==0;       
        figure; imagesc(imgtmp);
        h = imrect(gca,[110 10 400 400]);
        % Update position in title using newPositionCallback
        addNewPositionCallback(h,@(h) title(sprintf('(%1.0f,%1.0f)',h(1),h(2))));
        % Construct boundary constraint function
        fcn = makeConstrainToRectFcn('imrect',get(gca,'XLim'),get(gca,'YLim'));
        % Enforce boundary constraint function using setPositionConstraintFcn
        setPositionConstraintFcn(h,fcn);
        setColor(h,'r');
        %setString(h,'Plate center - then click')
        posctr = wait(h);
        mask= createMask(h);
        end
        close all %clean up 

        %% mask image for approval
        
        mask=double(mask);
        imgtmp=double(imgtmp);
        imgCP=(imgtmp.*mask);
        figure; imagesc(mask); colorbar

%if (allow_img) == 'Y' %Not optional- user needs to check 
    figure; 
        subplot (2,1,1); imagesc(imgCP);
        subplot (2,1,2); imagesc(imgtmp); hold on;  
   
%   cd(OutDir)
    %print('-dtiff', [ImageSetName, 'CrclMask'])
%end

%GoodCrop=input('is the crop good (y/n))', 's');
end
else %just  crop with the good parameters

        [xx,yy]=meshgrid(1:lng1,1:lng2); 
        rr=(xx-xctr).^2+(yy-yctr).^2;
        mask=rr<radSqr;
        mask=uint8(mask); %mask=double(mask);
        imgCP=(imgtmp.*mask);
        figure; imagesc(mask); colorbar
        
        if (allow_img) == 'Y'
        figure; 
        subplot (2,1,1); imagesc(imgCP);
        subplot (2,1,2); imagesc(imgtmp); hold on; 
        plot(xctr,yctr,'r*');
        theta=linspace(0,2*pi);
        plot(xctr+sqrt(radSqr)*cos(theta),yctr+sqrt(radSqr)*sin(theta), 'w');
        axis equal
        %cd(OutDir)
        %print('-dtiff', [ImageSetName, 'CrclMask'])
        end  
    
end

%center = [xctr,yctr];
%imwrite(imgCPtemp,[imgNum 'cropthresh.jpg']);
