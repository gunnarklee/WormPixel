%% PLOT FIGURES
    close all
    
    numFr =50; % frame count for the zoomed in view
    numFr2 =25;
    clims = [ -35 35 ];
    
    
    %full color full scale
    figure;imagesc(CurveMtx); colorbar
    saveas (gcf, [NameOut 'CurveMtxrfull'], 'pdf')
    
    MakeTricolMap
    %figure;imagesc(CurveMtx, clims)
    %colormap(RWBcMap2); colorbar
    %saveas (gcf, [NameOut 'CurveMtx'], 'pdf')
    
    figure;imagesc(CurveMtx(:,1:numFr), clims)
    colormap(RWBcMap2); colorbar
    saveas (gcf, [NameOut 'CurveMtx' num2str(numFr)], 'pdf')
    
    
    figure;imagesc(CurveMtx(:,1:numFr2), clims)
    colormap(RWBcMap2); colorbar
    saveas (gcf, [NameOut 'CurveMtx' num2str(numFr2)], 'pdf')
    
    figure;plot(centroid(:,1),centroid(:,2), '.k',...
        'LineWidth', 2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g', 'MarkerSize', 1);...
        title ('centroid position'); xlim ([1,size(img1,2)]); ylim([1, size(img1, 1)]);
    saveas (gcf, [NameOut 'PathTraveled'], 'pdf')
    
    %     figure;imagesc(CurveMtx)
    %     figure;plot(centroid(1:numFr,1),centroid(1:numFr,2), '--rs',...
    %         'LineWidth', 2, 'MarkerEdgeColor','k',...
    %         'MarkerFaceColor','g', 'MarkerSize', 1);...
    %         title ('centroid position'); xlim ([1,size(img1,2)]); ylim([1, size(img1, 1)]);
    %     saveas (gcf, [NameOut 'PathTraveled' num2str(numFr)], 'pdf')
    %
    %
    
    %figure;plot(centroid(:,1),centroid(:,2),'*b') ; title ('centroid position'); xlim ([1,450]); ylim([1,400]);
    
    %plot displacement
    figure;plot(time,distanceMv,'-r') ; title ('displacement vs. time');
    hold on
    axes('position', [.75 .15 .15 .15]);
    plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    saveas (gcf, [NameOut 'Distvstm'], 'pdf')
    
    %   figure;plot(time(1:numFr,:),distanceMv(1:numFr,:),'-r') ; title ('displacement vs. time');
    %   hold on
    %  axes('position', [.75 .15 .15 .15]);
    %   plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    %  saveas (gcf, [NameOut 'DistVstm' num2str(numFr)], 'pdf')
    
    %Cummulative Dist
    CumulDist=cumsum(distanceMv)
    figure;plot(time,CumulDist,'-r') ; title ('Cumldispl vs. time');
    hold on
    axes('position', [.75 .15 .15 .15]);
    plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    saveas (gcf, [NameOut 'CumlPathTraveled'], 'pdf')
    
    % figure;plot(time(1:numFr,:),CumulDist(1:numFr,:),'-r') ; title ('Cumldispl vs. time');
    % hold on
    % axes('position', [.75 .15 .15 .15]);
    % plot(varStruct.SpineData.SpineList(:,1), varStruct.SpineData.SpineList(:,2), 'r');
    %  saveas (gcf, [NameOut 'CumlPathTraveled' num2str(numFr)], 'pdf')
    
    %plot neck movement
    figure; plot(CurveMtx(3,1:numFr));
    saveas (gcf, [NameOut 'NeckMovement' num2str(numFr)], 'pdf')
    figure; plot(CurveMtx(3,:));
    saveas (gcf, [NameOut 'NeckMovement'], 'pdf')
    
    
    %     plot(SpineList{W}(:,1), SpineList{W}(:,2));
    %     hold on
    %     %figure; imagesc(CurveMtx);
    %
    
    %[wkSpnX,wkspnY]=ind2sub(size(SpineData.SpineList), find(SpineData.SpineList));
    %cmap=colormap(jet (size(wkSpnX,1))) %copper
    %WmImgPadcolor=mat2gray(WmImgPad)
    
    %WmImgPadcolor=(imoverlay (mat2gray(img1), SpineData.SpineList,  cmap(Pt,:)));
    
    %     figure; imshow(WmImgPadcolor, 'InitialMagnification', 400);
    %     hold on
    %     %addpoint to imgage with new color
    %     WmImgPadcolor=(imoverlay (mat2gray(WmImgPadcolor), skeleEND, cmap(Pt,:)));
end