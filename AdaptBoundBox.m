function [img1Masked] = AdaptBoundBox(boundingBox, PadPrc, MaskImage, img1, allow_img)
            switch MaskImage;
                case 'y' % LIMIT THE SEARCH AREA after the first image

                    BBmask = zeros(size(img1));
                    %Be sure that padded BB does not stretch the mask beyond the edges of img
                    %try % this fails somtimes
                        a=ceil(boundingBox(2)-boundingBox(4).*PadPrc);
                        b=ceil(boundingBox(2)+boundingBox(4).*2.*PadPrc);
                        c=ceil(boundingBox(1)-boundingBox(3).*PadPrc);
                        d=ceil(boundingBox(1)+boundingBox(3).*2.*PadPrc);
                        [a,b,c,d]=limit2Bounds(a,b,c,d,size(img1,1),size(img1,2));
                        BBmask(a:b,c:d)= ones;
                        [a,b,c,d, boundingBox];
                    %catch e1
                        %SpineData.FailPt= 'BBmask';
                    %    saveThis([ErrorDir filesep imageName(1:end-4), 'BBmaskError.mat'], varStruct)
                        %break
                    %end
                    % Be sure that padded BB does not stretch the mask beyond the edges of img
                    BBmask=uint8(BBmask(1:size(img1,1), 1:size(img1,2)));
                    img1Masked=(img1.*BBmask); %% pass the masked image into the particle analysis to avoit off target particles.
                    if (strcmpi (allow_img, 'y')) %display results
                        figure; imshow(uint8(img1Masked)); %PlotBoundBox(img1, boundingBox); % to check..
                    end
                case 'n'  %SECOND TRY, leave image UNMASKED
                    img1Masked=img1;
                    figure; imshow(uint8(img1Masked));
            end
end
    
