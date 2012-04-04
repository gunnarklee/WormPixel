 function Flipbook(path, imagelist)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for n=1:length(imagelist)

img1=imread([path filesep imagelist{n}]);
figure; imshow(img1);

end

