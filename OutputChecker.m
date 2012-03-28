function OutputChecker(finalfile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load(finalfile)
%load image and spine dat

%% SPECIFY PARAMETERS
OneWorm_CHR7Params

%Pad image in older versions I will need to reconsittute this padded image
if iscell(Imagesfilt)
    WmImg=Imagesfilt{1,1};
else
    WmImg=Imagesfilt;
end

%pad image
[WmImgPad] = padImg (WmImg, pad);


%%
%figure ;imshow(WmImgPad,'InitialMagnification', 400); title ('WmImgPad')
%figure; imshow(imoverlay (mat2gray(WmImgPad), SpineData.Pointlist(1),  [0, 0, 255]), 'InitialMagnification', 400);title ('skeleSH-shrunk');

%%

    MakeTricolMap

%%
%plot worm image image with spine and angles. 
    figure; subplot(2,3,[1:2.5,4:5.5]); imshow(WmImgPad, 'InitialMagnification', 400);
    hold on
    plot(SpineData.Pointlist(:,2), SpineData.Pointlist(:,1), 'g.', 'MarkerSize', 10)
    plot(SpineData.SpineList(:,2), SpineData.SpineList(:,1), 'b.', 'MarkerSize', 3)
    
    for pt=2:length(SpineData.Pointlist)-1
        text(SpineData.Pointlist(pt,2), SpineData.Pointlist(pt,1)-2,...
        num2str(SpineData.AngleLs(pt-1,1),2), 'Color',[0 1 0])
        pt
    end
  
%% Plot curverature column 
  
    clims=[-40 40]
    subplot(2,3,[3,6]); (imagesc(SpineData.AngleLs, clims)); colorbar; 
    colormap(RWBcMap)
    hold off
    %%

%>    cmap=Jet;
%>   for Pt=1:length(SpineData.SpineList)
%>   WmImgPadColor=(imoverlay (mat2gray(WmImgPad), SpineData.Pointlist(Pt,2), SpineData.Pointlist(Pt,1),  cmap(Pt,:)));
%>   end
%plot worm image image with spine and angles. 
%    figure; imshow(WmImgPadColor, 'InitialMagnification', 400);
 %   hold on
 %   plot(SpineData.Pointlist(:,2), SpineData.Pointlist(:,1), 'b*',
 %   'MarkerSize', 10)
 %RwB color map
%plot curverature column

end

