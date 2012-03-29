function [WmImgPad] = padImg (img, pad)
WmImgPad = zeros(size(img)+pad)
Rowpad=(pad/2+1:size(WmImgPad, 1)-pad/2)
ColPad=(pad/2+1:size(WmImgPad, 2)-pad/2)
WmImgPad(Rowpad, ColPad)=img
%figure ;imshow(WmImgPad)
end