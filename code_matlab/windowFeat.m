function [EpochFeat,StEpochTime,EpochData]=windowFeat(data,fs,opt)
%opt.twin  Window on which Breathing features are estimated
%opt.twinMove Window for  window slide 
%feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
%          skew_mean kurt_mean entro...
%          per_power  cycle fmax calibDel];
%  EpochFeat = num of epoch * feature num * channel
%  StEpochTime = start time of each epoch (sec)  based on Start time(full time)
%  Data first normalized to [-1,1] in each epoch


twinMove=opt.twinMove;
twin=opt.twin;



StEpochTime=0:twinMove:size(data,1)/fs-twin;  %90s epoch number   
StEpochTime=StEpochTime';

EpochData=zeros(length(StEpochTime),twin*fs,size(data,2));% epoch num * epoch samplePoint * channel
EpochFeat=zeros(length(StEpochTime),opt.featNum,size(data,2));
for i=1:length(StEpochTime)
     idx=[1+StEpochTime(i)*fs:(StEpochTime(i)+twin)*fs];
     EpochData(i,:,:)=data(idx,:);
     EpochDataTemp=data(idx,:);
     
     EpochFeat(i,:,:)=FeatInEpoch(EpochDataTemp,fs,opt);
     
    
end
end