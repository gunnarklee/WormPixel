function [vertLn,horzLn] = plot4Stacktiff (plotpoints, CentrComp, plotcol, I, F, textls, PlotMd, LayerNm, boundingBox, BBratio)


Smoothstart=50
%plot4Stacktiff, L. Parsons and G. kleemnann 7_15_11
%function to automatically make a stacked image tiff uisng imshow and/or imagesc 
%it plots title, subtitle and socred points options

%INPUT VALUES
% plotpoints,-y/n
% plotcol, - columns in CentrComp to plot [3,4]
% I, 
% F, 
% textls - CurImgNm; CurImgTitle; Assaytime; Measur
% CentrComp,   - point matrox referred to by plotcol 
% PlotMd, 
% LayerNm

%EXAMPLE
%Exisiting matricies

%Img_Propfilt %filtered particles with location in col 6,7
%img1   %image in matirx
%img2

%Nm = 'time t';
%Ti = 'Image 1, NO POINTS'
%F = figure('Visible', 'on');  % Create a new figure without displaying it
%plotP = 'n';
%plot4Stacktiff (plotP, 6,7, uint8(img1), F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
%I = imread('Image.tif');
%imwrite(I, 'stack.tif');

%Nm = 'time t+1';;
%Ti = 'Image 2, NO POINTS';
%F = figure('Visible', 'on');  % Create a new figure without displaying it
%plotP = 'n';
%plot4Stacktiff (plotP, 6,7, uint8(img2), F, Nm, Ti, Img_Propfilt,'imshow', 'Image.tif');
%I = imread('Image.tif');
%imwrite(I, 'stack.tif', 'writemode', 'append');


I_sz = size(I); % Get the size of the image
% Factor to use if image is too large for screen 
% NOTE: Using this may require that plotted elements size be adjusted
N = 1; 
% Strip excess edges (only enough space to show image)
set(F, 'Units', 'pixels', 'Position', [32, 32, I_sz(2) / N, I_sz(1) / N]); 
set(gca, 'Units', 'normalized', 'Position', [0,0,1,1]); 
set(gca, 'Units', 'points');

% Display image and plot on top of it (replace this with your code)
%imagesc(I)

    if strcmpi('imshow', PlotMd)
    imshow(I) %
    else
    imagesc(I)
    end
    
hold on
xlims=get(gca,'xlim');
ylims=get(gca,'ylim');
x=xlims(1):1:xlims(2);

    %if strcmpi('y', plotpoints)
    %end

%make cross hairs for bounding box
horzLn=[boundingBox(1):boundingBox(1)+boundingBox(3)]';
PosH=repmat(CentrComp(size(CentrComp, 1),2), 1, size(horzLn,1))';
horzLn=[PosH,horzLn];
vertLn=[boundingBox(2):boundingBox(2)+boundingBox(4)]';
PosV=repmat(CentrComp(size(CentrComp, 1),1), 1, size(vertLn,1))';
vertLn=[PosV, vertLn];

if size(CentrComp, 1) < Smoothstart
plot (CentrComp(:,plotcol(1)), CentrComp(:, plotcol(2)),'b*', CentrComp(:,plotcol(1)), CentrComp(:, plotcol(2)),'-g')
else
plot (CentrComp(:,plotcol(1)), CentrComp(:, plotcol(2)),'b*', CentrComp(:,plotcol(3)), CentrComp(:, plotcol(4)),'-r')   
end

plot (horzLn(:,2),horzLn(:,1),'-g','LineWidth',2)
    plot (vertLn(:,1),vertLn(:,2),'-g','LineWidth',2)

CurImgTitle=textls{1} %text to be printed
CurImgNm=textls{2}
Assaytime=textls{3}
Measur=textls{4}

text(50,10, ['\fontsize{25}', '\color{magenta}',CurImgTitle]); 
text(50,450, ['\fontsize{16}', '\color{magenta}',CurImgNm]) ;
text(50,425, ['\fontsize{16}', '\color{magenta}', Assaytime]);
text(50,400, ['\fontsize{16}', '\color{magenta}', Measur]);
          
    %>>     axes ('position', [.65, .1, .3, .3]);
    %>>     plot(BBratio(2,:),BBratio(1,:));
    %>>    title('Major/minor axis Vs image number');


% Poisition figure for "printing" with no border
set(F, 'Units', 'points', 'PaperUnits', 'points', 'PaperPositionMode', 'auto');
% Print (save) figure as bitmapped image, resizing by factor N
screen_DPI = get(0, 'ScreenPixelsPerInch');
print(F, '-dtiff', sprintf('-r%d', N * screen_DPI), LayerNm);
%print(F, '-dpng', sprintf('-r%d', N * screen_DPI), 'pout_with_plot.png');
% Close figure
%sclose(F)
