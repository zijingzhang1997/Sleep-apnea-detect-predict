function [EpochFeat]=windowFeat_HF(data,fs,opt)
%opt.twin  Window on which Breathing features are estimated
%opt.twinMove Window for  window slide 
%High frequency feature=[psd(1) psd(2) psd(3) psd(4)]
      


twinMove=opt.twinMove;
twin=opt.twin;



StEpochTime=0:twinMove:size(data,1)/fs-twin;  %90s epoch number   
StEpochTime=StEpochTime';

EpochData=zeros(length(StEpochTime),twin*fs,size(data,2));% epoch num * epoch samplePoint * channel
EpochFeat=zeros(length(StEpochTime),opt.featNumHF,size(data,2));
for i=1:length(StEpochTime)
     idx=[1+StEpochTime(i)*fs:(StEpochTime(i)+twin)*fs];
     EpochData(i,:,:)=data(idx,:);
     EpochDataTemp=data(idx,:);
     
     EpochFeat(i,:,:)=FeatInEpoch_HF(EpochDataTemp,fs,opt);
     
    
end
end

function featAll=FeatInEpoch_HF(EpochData,fs,opt)
featAll=zeros(opt.featNumHF,size(EpochData,2)); 

for i=1:size(EpochData,2)
    EpochData_row=EpochData(:,i);
    
    featAll(:,i)=featExtract_HF(EpochData_row,fs,opt);
    
end  

end

    

function feature=featExtract_HF(Data,fs,opt)    
    





    %% Find features in each cycle 


[pxx,f] = periodogram(Data,hamming(length(Data)),length(Data),fs,'power');
f(1:5)=[];pxx(1:5)=[];

pTot = bandpower(pxx,f,'psd');
opt.HFband(1,1)=f(1);
pBand(1) = bandpower(pxx,f,[opt.HFband(1,1) opt.HFband(1,2)],'psd');
pBand(2) = bandpower(pxx,f,[opt.HFband(2,1) opt.HFband(2,2)],'psd');
pBand(3) = bandpower(pxx,f,[opt.HFband(3,1) opt.HFband(3,2)],'psd');
pBand(4) = bandpower(pxx,f,[opt.HFband(4,1) opt.HFband(4,2)],'psd');




pBand(:) = pBand(:).*(100/pTot);



% find psd in specific bandwidth
f_b1=f(find(f>opt.HFband(1,1)&f<opt.HFband(1,2)));
pxx_b1=pxx(find(f>opt.HFband(1,1)&f<opt.HFband(1,2)));
fmax(1)=f_b1(find(pxx_b1==max(pxx_b1)));

f_b2=f(find(f>opt.HFband(2,1)&f<opt.HFband(2,2)));
pxx_b2=pxx(find(f>opt.HFband(2,1)&f<opt.HFband(2,2)));
fmax(2)=f_b2(find(pxx_b2==max(pxx_b2)));

f_b3=f(find(f>opt.HFband(3,1)&f<opt.HFband(3,2)));
pxx_b3=pxx(find(f>opt.HFband(3,1)&f<opt.HFband(3,2)));
fmax(3)=f_b3(find(pxx_b3==max(pxx_b3)));

f_b4=f(find(f>opt.HFband(4,1)&f<opt.HFband(4,2)));
pxx_b4=pxx(find(f>opt.HFband(4,1)&f<opt.HFband(4,2)));
fmax(4)=f_b4(find(pxx_b4==max(pxx_b4)));



     
feature=[pBand fmax];
     
feature(isnan(feature))=0;




    
end