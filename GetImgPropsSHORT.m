function [imgBWL, F, Image_PropertiesAll] = GetImgProps (BW, allow_img)
%% GetImgProps
%G.kleemann 2/6/08
%Input base thresholded image > BW=BWtp
%output imgBWL=imgBWLtp, Image_PropertiesAll
%% use BW label to get objects

%BW=imcomplement(BW);
%[AlabeldAll]=bwconncomp(BW,4); %faster
[AlabeldAll,gnumAll]=bwlabel(BW,4);
imgBWL=label2rgb(AlabeldAll, @jet, 'k', 'shuffle');

%>close all
%% calculate and record properties and 

%>>>F=regionprops(AlabeldAll,'Area','Eccentricity','MajorAxisLength','MinorAxisLength',...
%>>>    'Orientation', 'centroid','BoundingBox', 'Extent','EquivDiameter','EulerNumber',...
%>>>    'Perimeter','PixelList', 'Solidity', 'FilledImage', 'image');
 


F=regionprops(AlabeldAll,'Area','MajorAxisLength','MinorAxisLength','centroid','BoundingBox', 'image');


 %F=regionprops(AlabeldAll, 'all');

 %REMOVING THIS LOOP WILL SPEED THINGS UP!
 %for object_num=1:gnumAll
  %      [r,c]=find(AlabeldAll==object_num);
   %     Image_PropertiesAll(object_num,1:4)=[mean(r),mean(c),length(r),object_num,]; %[Xcordinate Ycoordinate Worm_or_dust_Size(length of vector)]
 %end

%% 
 
%cd(Programs)
%>ip=input_points;
%>bp=base_points;

%Make matrix of eccentricity and certroid values (EccMtx), this should be concatenated to Image_properties for filtering!
%fortunately both lists are in the object numnber order, so they can be
%concatinated   - compare centroid columns 1,2 to 6,7 
Ecc = cat(1, []);
Ctrd = cat(1, F.Centroid);
Area= cat(1, F.Area);
MajAx= cat(1, F.MajorAxisLength);
MinAx= cat(1, F.MinorAxisLength);
%MajtoMinAx= cat(1, (F.MajorAxisLength/F.MinorAxisLength)); % this may be a worm typical #, but curled worm will throw this off
BB= cat(1, F.BoundingBox);
Ext= cat(1, []); % of bounding box filled, perhaps combine with mintomaxratio
EqD= cat(1, []);
EuN= cat(1, []);
Perim= cat(1, []);
objnum=[1:gnumAll]';
perim_area=[];%Perim./Area;
majAx_minAx=MajAx./MinAx;
majAx_minAx_ext=[];%(MajAx./MinAx)./Ext;

Image_PropertiesAll = cat(2, Ctrd, Area, objnum, Ecc, Ctrd, MajAx, MinAx, Area, BB, Ext, EqD, EuN, Perim, perim_area, majAx_minAx, majAx_minAx_ext);
%first two items (three columns) are to retain spacing

if (strcmpi (allow_img, 'y'))
    figure;    
        %subplot (3,3,1); hist (Image_PropertiesAll (:,5)); title ('Eccentr 1=straight')
        subplot (3,3,2); hist (Image_PropertiesAll (:,8)); title ('MajAxis')
        subplot (3,3,3); hist (Image_PropertiesAll (:,9)); title ('MinAxis')
        subplot (3,3,4); hist (Image_PropertiesAll (:,10)); title ('Area')
       % subplot (3,3,5); hist (Image_PropertiesAll (:,15)); title ('Extent')
        %subplot (3,3,6); hist (Image_PropertiesAll (:,16)); title ('equivalentDiameter')
        %subplot (3,3,7); hist (Image_PropertiesAll (:,17)); title ('EulerNumber')
        %subplot (3,3,8); hist (Image_PropertiesAll (:,18)); title ('Perimiter/area')
        subplot (3,3,9); hist (Image_PropertiesAll (:,20)); title ('maj_min_extnt')
end    