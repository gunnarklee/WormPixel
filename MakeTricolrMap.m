%% making tricolor maps    
    RWBcMap = [1 0 0; 1 1 1 ; 0 0 1] 
    RWBcMap2=zeros(64,3);
    %RWBcMap2(1:32,1) = sort([0:1/31:1]', 'descend')
    %RWBcMap2(33:64,3) = sort([0:1/31:1]', 'ascend')
    
    RWBcMap2(1:32,1) = ones
    RWBcMap2(33:64,3) = ones
    
    RWBcMap2(1:32,2) = sort([0:1/31:1]', 'ascend')
    RWBcMap2(33:64,2) = sort([0:1/31:1]', 'descend')
    
      figure; subplot (1,2,1); imagesc(RWBcMap);...
              subplot (1,2,2); imagesc(RWBcMap2);