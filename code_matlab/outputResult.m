function  [feature,label]=outputResult(flag,EpochFeat,EpochFeatPSGSpO2,labelEpoch,opt,validEpochIdx)

idx0=find(validEpochIdx==0 );
featNum=size(EpochFeat,2);
featNum2=size(EpochFeatPSGSpO2,2);

if opt.ncs ==1  %NCS need to select true epoch
    flag(idx0,:)=0;
    [r,c]=find(flag==1);
    
    for i=1:length(r)   % total number of true NCS epoch
        feature(i,1:featNum)=EpochFeat(r(i),:,c(i));
        feature(i,featNum+1:featNum+featNum2)=EpochFeatPSGSpO2(r(i),:);
        label(i)=labelEpoch(r(i));
        
    
    end
    
    
end
if opt.ncs ==2  % PSG choose the same unique epoches and use the optimal channel
    flag(idx0,:)=0;
    [r,~]=find(flag==1);
    r=unique(r);
    psgIdx=opt.psgIdx;
    
    for i=1:length(r)
        feature(i,1:featNum)=EpochFeat(r(i),:,psgIdx(r(i)));
        feature(i,featNum+1:featNum+featNum2)=EpochFeatPSGSpO2(r(i),:);
        label(i)=labelEpoch(r(i));
        
    end
    
    
end
if opt.ncs==0 %PSG  no need to select true epoch . use all epoch
    num=1;
    psgIdx=opt.psgIdx;
    for i=1:size(EpochFeat,1)
        
        
        feature(num,1:featNum)= EpochFeat(i,:,psgIdx(i));
        feature(num,featNum+1:featNum+featNum2)=EpochFeatPSGSpO2(i,:);
        label(num)=labelEpoch(i);
        num=num+1;
        
    end
    feature(idx0,:)=[];
    label(idx0)=[];
end

label=label';

end