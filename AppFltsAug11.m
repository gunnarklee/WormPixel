%%AppFlts


%set up filters to identify worms - logical matricies of ones and zeros...
szFltL= double(Image_PropertiesAll(:,10)>LowLim); %make size filter 
szFltU= double(Image_PropertiesAll (:,10)<UpLim);
EccFltL= double(Image_PropertiesAll (:,5)>eccL); %rows less than 1280 are in chanel 1
EccFltU= double(Image_PropertiesAll (:,5)<eccU);% the eccentricity of the spot should be less than .8 (usually ~.5)
MajAxFltL= double(Image_PropertiesAll (:,8)>MajAxL); %rows less than 1280 are in chanel 1
MajAxFltU= double(Image_PropertiesAll (:,8)<MajAxU);% the eccentricity of the spot should be less than .8 (usually ~.5)
MinAxFltL= double(Image_PropertiesAll (:,9)>MinAxL); %rows less than 1280 are in chanel 1
MinAxFltU= double(Image_PropertiesAll (:,9)<MinAxU);% the eccentricity of the spot should be less than .8 (usually ~.5)
ExtFltL= double(Image_PropertiesAll (:,15)>ExtL); %rows less than 1280 are in chanel 1
ExtFltU= double(Image_PropertiesAll (:,15)<ExtU);   
TotAxFltU = (Image_PropertiesAll (:,8))+(Image_PropertiesAll (:,9))<TotAxU;
TotAxFltL = (Image_PropertiesAll (:,8))+(Image_PropertiesAll (:,9))>TotAxL;

