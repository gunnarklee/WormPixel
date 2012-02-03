function [output_args ] = saveImageToStack(image, filename, varargin)
%SAVEIMAGETOSTACK Summary of this function goes here
%   Detailed explanation goes here

% Parse arguments
p = inputParser;
p.addRequired('image', @isnumeric);
p.addRequired('filename', @ischar);
p.addParamValue('title', 'Title', @ischar);
p.addParamValue('image_name', 'Image Name', @ischar);
p.addParamValue('plot_points', 0, @isnumeric);
p.addParamValue('scale_image', false, @islogical);
p.addParamValue('display_image', 'off', @ischar);
p.addParamValue('CentrComp', [1,1,1,1],  @isnumeric)
p.addParamValue('boundingBox', [1,1,1,1] , @isvector)
p.addParamValue('plotcol', [1,2,3,4], @isvector)

p.parse(image, filename, varargin{:})

% Create figure
F = figure('Visible', p.Results.display_image);  % Create a new figure without displaying it

% Get the size of the image
I_sz = size(image);

% Factor to use if image is too large for screen
% NOTE: Using this may require that plotted elements size be adjusted
N = 1;

% Strip excess edges (only enough space to show image)
set(F, 'Units', 'pixels', 'Position', [32, 32, I_sz(2) / N, I_sz(1) / N]);
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1]);
set(gca, 'Units', 'points');

% Display image
if p.Results.scale_image
    imagesc(image)
else
    imshow(image)
end


%% Plot data
%make cross hairs for bounding box
Smoothstart=50
horzLn=[p.Results.boundingBox(1):p.Results.boundingBox(1)+p.Results.boundingBox(3)]';
PosH=repmat(p.Results.CentrComp(size(p.Results.CentrComp, 1),2), 1, size(horzLn,1))';
horzLn=[PosH,horzLn];
vertLn=[p.Results.boundingBox(2):p.Results.boundingBox(2)+p.Results.boundingBox(4)]';
PosV=repmat(p.Results.CentrComp(size(p.Results.CentrComp, 1),1), 1, size(vertLn,1))';
vertLn=[PosV, vertLn];

hold on 
if size(p.Results.CentrComp, 1) < Smoothstart
    plot (p.Results.CentrComp(:,p.Results.plotcol(1)),...
        p.Results.CentrComp(:, p.Results.plotcol(2)),'b*',...
        p.Results.CentrComp(:,p.Results.plotcol(1)),...
        p.Results.CentrComp(:, p.Results.plotcol(2)),'-g')
else
    plot (p.Results.CentrComp(:,p.Results.plotcol(1)),...
        p.Results.CentrComp(:, p.Results.plotcol(2)),'b*',...
        p.Results.CentrComp(:,p.Results.plotcol(3)),...
        p.Results.CentrComp(:, p.Results.plotcol(4)),'-r')
end

plot (horzLn(:,2),horzLn(:,1),'-g','LineWidth',2)
plot (vertLn(:,1),vertLn(:,2),'-g','LineWidth',2)


%text to be printed
% CurImgTitle=textls{1}
% CurImgNm=textls{2}
% Assaytime=textls{3}
% Measur=textls{4}
%
% text(50,10, ['\fontsize{25}', '\color{magenta}',CurImgTitle]);
% text(50,450, ['\fontsize{16}', '\color{magenta}',CurImgNm]) ;
% text(50,425, ['\fontsize{16}', '\color{magenta}', Assaytime]);
% text(50,400, ['\fontsize{16}', '\color{magenta}', Measur]);


% ****Plot points on top of image > not used yet 
if p.Results.plot_points
    plot (p.Results.plot_points(:,1), p.Results.plot_points(:,2),'r*', 'MarkerSize', 3)
end
text(50,50, ['\fontsize{25}', '\color{magenta}',p.Results.title])
text(50,450, ['\fontsize{16}', '\color{magenta}',p.Results.image_name])

% Poisition figure for "printing" with no border
set(F, 'Units', 'points', 'PaperUnits', 'points', 'PaperPositionMode', 'auto');

% Print (save) figure as bitmapped image, resizing by factor N
screen_DPI = get(0, 'ScreenPixelsPerInch');
print(F, '-dpng', sprintf('-r%d', N * screen_DPI), [filename, '.tmp']);

% Close figure
if (strcmpi(p.Results.display_image, 'n'))
    close(F)
end

% Append to file
I = imread([filename, '.tmp']);
delete([filename, '.tmp']);
[IND,map] = rgb2ind(I,256);
if exist(filename, 'file')
    imwrite(IND, map, filename, 'gif', 'writemode', 'append');
else
    imwrite(IND, map, filename, 'gif');
end

end