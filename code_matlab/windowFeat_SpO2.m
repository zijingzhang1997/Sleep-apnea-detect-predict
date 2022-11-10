function [EpochFeat,StEpochTime,EpochData]=windowFeat_SpO2(data,fs,opt)
%opt.twin  Window on which Breathing features are estimated
%opt.twinMove Window for  window slide 
%feature=[mean(SpO2) std(SpO2) percent(SpO2<threshold) minimum (SpO2)]
      


twinMove=opt.twinMove;
twin=opt.twin;
opt.lowTh;


StEpochTime=0:twinMove:size(data,1)/fs-twin;  %90s epoch number   
StEpochTime=StEpochTime';

EpochData=zeros(length(StEpochTime),twin*fs,size(data,2));% epoch num * epoch samplePoint * channel
EpochFeat=zeros(length(StEpochTime),4,size(data,2));
for i=1:length(StEpochTime)
     idx=[1+StEpochTime(i)*fs:(StEpochTime(i)+twin)*fs];
     EpochData(i,:)=data(idx);
     EpochDataTemp=data(idx);
     
     EpochFeat(i,1)=mean(EpochDataTemp);
     EpochFeat(i,2)=std(EpochDataTemp);
     EpochFeat(i,3)=min(EpochDataTemp);
     EpochFeat(i,4)=length(find(EpochDataTemp<=opt.lowTh*100))/(length(EpochDataTemp));
     
end

end