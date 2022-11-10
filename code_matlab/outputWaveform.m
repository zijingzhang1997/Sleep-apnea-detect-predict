% Output waveform data directly, instead of features 
% for CNN training purpose 
% feature size = case  * waveform sample  * channel 


%opt.ncs =1 :NCS need to select true epoch 
% =2 PSG2 choose the same unique epoches and use the optimal channel 
% =0 PSG  no need to select true epoch . use all epoch


% 2 channels 
% ch1: NCS/PSG waveform  do normalize [-1,1]? 
% ch2: SpO2 waveform   do normalize *0.01? 

function  [feature,label]=outputWaveform(flag,EpochData,EpochDataPSGSpO2,labelEpoch,opt,validEpochIdx)

idx0=find(validEpochIdx==0 );
featNum=size(EpochData,2);


if opt.ncs ==1  %NCS need to select true epoch
    flag(idx0,:)=0;
    [r,c]=find(flag==1);
    
    for i=1:length(r)
        tempData1=EpochData(r(i),:,c(i));
        
        tempData2=EpochDataPSGSpO2(r(i),:);
        tempData1 = rescale(tempData1,-1,1);
        tempData2 = tempData2.*0.01;
        feature(i,1:featNum,1)=tempData1;
        feature(i,1:featNum,2)=tempData2;
        
        label(i)=labelEpoch(r(i));
        
        
        
    end
    
    
end
if opt.ncs ==2  % PSG choose the same unique epoches and use the optimal channel
    flag(idx0,:)=0;
    [r,~]=find(flag==1);
    r=unique(r);
    psgIdx=opt.psgIdx;
    
    for i=1:length(r)
  
        tempData1=EpochData(r(i),:,psgIdx(r(i)));
        tempData2=EpochDataPSGSpO2(r(i),:);
        tempData1 = rescale(tempData1,-1,1);
        tempData2 = tempData2.*0.01;
        feature(i,1:featNum,1)=tempData1;
        feature(i,1:featNum,2)=tempData2;
        
        label(i)=labelEpoch(r(i));
    end
    
    
end
if opt.ncs==0 %PSG  no need to select true epoch . use all epoch
    num=1;
    psgIdx=opt.psgIdx;
    for i=1:size(EpochData,1)
        
        tempData1=EpochData(i,:,psgIdx(i));
        tempData2=EpochDataPSGSpO2(i,:);
        tempData1 = rescale(tempData1,-1,1);
        tempData2 = tempData2.*0.01;
        feature(num,1:featNum,1)=tempData1;
        feature(num,1:featNum,2)=tempData2;
       
        
        label(num)=labelEpoch(i);
        num=num+1;
        
    end
    feature(idx0,:,:)=[];
    label(idx0)=[];
    
end

label=label';

end