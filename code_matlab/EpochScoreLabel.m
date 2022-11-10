function [labelEpoch,labelProp,validEpochIdx] =EpochScoreLabel(score_all,score,StEpochTime,opt,fs)

% rate=[apneaRaio(in epoch)  apeanRatio(for all seconds)];
% class 1. label 0 = score 0-4 (no apnea)   label 1 =  score 5-10 (apnea)
% threshold=opt.labelratio

% score   col 1: time sample      col 2:end : label 0-12  (0/1)
%  validEpochIdx (0/1)  invalid means score 0 > 0.5*length
validEpochIdx=ones(length(StEpochTime),1);
labelEpoch=zeros(length(StEpochTime),1);
ratioTh=opt.labelRaio;
ratioThSnore=opt.ratioThSnore;
    
if opt.classMethod==1
 for i=1:length(StEpochTime)
    idx=[StEpochTime(i)*fs+1:(StEpochTime(i)+opt.twin)*fs];

    EpochTemp=score_all(idx,2:end);
    for j=1:12   %label number 0-12    length (num)=13; label=numIndex -1;
        num(j)=length(find(EpochTemp(:,j)==1));
    end
    numLabel(1)=num(5);    %  LABEL 1 :snore  score=5
    numLabel(2)=sum(num(6:9));   %LABEL 2 :RERA score=7 + arousal score=6 hyponea  score=9 + O2 desaturation score=8
   
    numLabel(3)=num(10:11);   %LABEL 3 :OSA  score=10 
  
    numLabel(4)=num(12);   %LABEL 4 :CSA  score=12 
    ratioTh=opt.labelRaio;
    
    for m=1:length(numLabel)
        if numLabel(m)>=ratioTh*length(EpochTemp)
            labelEpoch(i)=m;
        end
    end
    
    
    
    if sum(num)<= 0.5 * length(EpochTemp)
        validEpochIdx (i) =0;
    end
  end         

end
%% 

if opt.classMethod==2
         
  for i=1:length(StEpochTime)
    idx=[StEpochTime(i)*fs+1:(StEpochTime(i)+opt.twin)*fs];

    EpochTemp=score_all(idx,2:end);
    EpochTemp2=score(idx,2);
    EpochTemp3=max(EpochTemp,2);
    for j=1:12   %label number 0-12    length (num)=13; label=numIndex -1;
        num(j)=length(find(EpochTemp(:,j)==1));
    end
    numLabel(1)=num(5);    %  LABEL 1 :snore  score=5
    numLabel(2)=sum(num(6:7));   %LABEL 2 :RERA score=7 + arousal score=6
    
    numLabel(3)=sum(num(8:9));   %LABEL 3 :hyponea  score=9 + O2 desaturation score=8
  

    numLabel(4)=num(10);   %LABEL 4 :OSA  score=10 
    numLabel(5)=num(11);   %LABEL 5 :MA  score=11 
    numLabel(6)=num(12);   %LABEL 6 :CSA  score=12 
    
    for m=2:length(numLabel)
        if numLabel(m)>=ratioTh*length(EpochTemp)
            labelEpoch(i)=m;
        end
    end
    
    if numLabel(1)>=ratioThSnore*length(EpochTemp)
            labelEpoch(i)=1;
    end
    
    
    
    if sum(num)<= 0.5 * length(EpochTemp)
        validEpochIdx (i) =0;
    end
  end
end
%% 


numapnea=length(find(labelEpoch~=0));
apneaRaio=numapnea/length(StEpochTime);

[C,ia,ic] = unique(labelEpoch);
a_counts = accumarray(ic,1);
label_counts_epoch = [C, a_counts];

[C,ia,ic] = unique(score(:,2));
a_counts = accumarray(ic,1);
label_counts_all = [C, a_counts];



labelProp.apneaRaioEpoch=apneaRaio;
labelProp.apneaRaioAll=length(find(score_all(:,6:end)==1))/length(score_all);
labelProp.label_counts_epoch=label_counts_epoch;
labelProp.label_counts_all=label_counts_all;
end