function [HeadPos] = GetPoint (img, startpos, mssg)
% INPUT = img, startpos
%   img = Image used to position crop ellipse
%   startpos = start position of the crop ellipse [x, y, w, h]
% OUTPUT = [mask, posctr]
%   mask = binary matrix where 0 outside of selected circle, 1 inside
%   posctr = x,y coordinates of image center

%CropCircleAdptAuto3 Gunnar Kleemann (6/19/08)
%version 3 adds optional image supression
%Circular masking to select 4 regions

disp(mssg);
done = 0;
while ~done;
    
    f = figure; imagesc(img);
    h = impoint(gca,startpos);
    
    % Update position in title using newPositionCallback
    addNewPositionCallback(h,@(h) title(sprintf('(%1.0f,%1.0f)',h(1),h(2))));
    % Construct boundary constraint function
    fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
    % Enforce boundary constraint function using setPositionConstraintFcn
    setPositionConstraintFcn(h,fcn);
    setColor(h,'r');
    %setString(h,'Plate center - then click')
    done = wait(h);
    %mask= createMask(h);
    HeadPos = getPosition(h);
end
close(f) %clean up
