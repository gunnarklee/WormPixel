%PrintObjParams.m
%INPUT (Image_PropertiesList, ctX, ctY, TxtCol))

plot (Image_PropertiesList(:,ctX), Image_PropertiesList(:,ctY),'rO', 'MarkerSize', 1)

for rowN = 1:length(Image_PropertiesList(:,1))
    text(Image_PropertiesList(rowN,ctX),...
    Image_PropertiesList(rowN,ctY),...
    num2str(Image_PropertiesList(rowN,TxtCol)),'Color', NumCol, 'FontSize', FntSz)
end