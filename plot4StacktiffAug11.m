function plot4Stacktiff (plotpoints, xcol, ycol, I, F, CurImgNm, CurImgTitle, pnts, PlotMd, LayerNm)

%plot4Stacktiff, L. Parsons and G. kleemnann 7_15_11
%function to automatically make a stacked image tiff uisng imshow and/or imagesc 
%it plots title, subtitle and socred points options

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
hold on

    if strcmpi('y', plotpoints)
    plot (pnts(:,xcol), pnts(:,ycol),'r*', 'MarkerSize', 3)
    end
text(50,50, ['\fontsize{25}', '\color{magenta}',CurImgTitle]) 
text(50,450, ['\fontsize{16}', '\color{magenta}',CurImgNm]) 

% Poisition figure for "printing" with no border
set(F, 'Units', 'points', 'PaperUnits', 'points', 'PaperPositionMode', 'auto');

% Print (save) figure as bitmapped image, resizing by factor N
screen_DPI = get(0, 'ScreenPixelsPerInch');
print(F, '-dtiff', sprintf('-r%d', N * screen_DPI), LayerNm);
%print(F, '-dpng', sprintf('-r%d', N * screen_DPI), 'pout_with_plot.png');
% Close figure
%sclose(F)
