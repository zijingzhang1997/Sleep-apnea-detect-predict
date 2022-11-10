function  []=EpochComp(PSGdata,NCSdata,fs,StEpochTime,opt,)
StEpochTime=StEpochTimeNCS;
twin=opt.twin;
cor11=zeros(length(StEpochTime),1);


for i=1:length(StEpochTime)
     idx=[1+StEpochTime(i)*fs:(StEpochTime(i)+twin)*fs];
     
     EpochTempNCS=NCSdata(idx,:);
     EpochTempPSG=PSGdata(idx,:);
     %rescale to [0,1] in every epoch
     EpochTempNCS = rescale(EpochTempNCS,'InputMin',min(EpochTempNCS,[],1),'InputMax',max(EpochTempNCS,[],1));
     EpochTempPSG = rescale(EpochTempPSG,'InputMin',min(EpochTempPSG,[],1),'InputMax',max(EpochTempPSG,[],1));
     [c,~] = xcorr(EpochTempPSG(:,1),EpochTempNCS(:,1),fs*10,'normalized');
     cor11(i)=max(c);
     
end


end