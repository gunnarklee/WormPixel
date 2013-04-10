function [row, mindiff]=CloseCentr(pos, Image_PropertiesAll)

distancediff=[];
for Ptnm=1:size(Image_PropertiesAll, 1)
    
    %get each distance between the seleced position and each centroid
    currdist=pdist2(Image_PropertiesAll(Ptnm,1:2),pos);
    [row,col]=find(distancediff==min(distancediff));
    distancediff=[distancediff;currdist];
end
%row and dist for the closest particle
[row,col]=find(distancediff==min(distancediff));
mindiff=distancediff(row,:);
end