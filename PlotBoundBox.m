function PlotBoundBox(img, boundingBox)

horzLn=[boundingBox(1):boundingBox(1)+boundingBox(3)]';
PosH=repmat(boundingBox(2),1, size(horzLn,1))';
PosH2=repmat(boundingBox(2)+boundingBox(4),1,size(horzLn,1))';
horzLn1=[PosH,horzLn];
horzLn2=[PosH2,horzLn];   

VertLn=[boundingBox(2):boundingBox(2)+boundingBox(4)]';
PosV=repmat(boundingBox(1),1, size(VertLn,1))';
PosV2=repmat(boundingBox(1)+boundingBox(3),1,size(VertLn,1))';
vertLn1=[VertLn, PosV];
vertLn2=[VertLn, PosV2];

figure; imshow(uint8(img))
hold on
plot (horzLn1(:,2),horzLn1(:,1),'-g','LineWidth',2)
plot (horzLn2(:,2),horzLn2(:,1),'-g','LineWidth',2)
plot (vertLn1(:,2),vertLn1(:,1),'-g','LineWidth',2)
plot (vertLn2(:,2),vertLn2(:,1),'-g','LineWidth',2)

end