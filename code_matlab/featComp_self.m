function [flag,property,optIdx] =featComp_self(EpochFeat1,EpochFeat2,...
    opt1,opt2,opt3,StEpochTime,fs,labelEpoch)
%feat 1 reference   feat 2 compare object
%feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)...
%          skew_mean kurt_mean entro...
%          per_power fmax];
%  compare [mean(br)1 mean(in)5 mean(ex)7  entro per_power12 ]

%compare differently accroding to apnea/ noapnea
flag=zeros(size(EpochFeat2,1,3));
MeetLimit=[0 0 0];NoMeetDiff=[0 0 0];
falseNum=[0 0 0];trueNum=[0 0 0]; NCSch=[0 0 0 0];NoMeetP=[0 0 0 0];
labelMark=[0 0 0];

%select the optimal chennel of PSG by least std. Record the idx
for n=1:size(EpochFeat1,1)
    for m=1:size(EpochFeat1,3)
        stdTemp(m)=(EpochFeat1(n,2,m)/EpochFeat1(n,1,m)+EpochFeat1(n,4,m)/EpochFeat1(n,3,m)...
            +EpochFeat1(n,6,m)/EpochFeat1(n,5,m)+EpochFeat1(n,8,m)/EpochFeat1(n,7,m));
    end
    
    [~,ind]=min(stdTemp);
    EpochFeat1_mean(n,:)=EpochFeat1(n,:,ind);
    optIdx(n)=ind;
end

for m=1:size(EpochFeat2,3)
    EpochFeat2_diff(:,:,m)=EpochFeat1_mean-EpochFeat2(:,:,m);
    
end




for j=1:size(EpochFeat2,3)
    for i=1:length(EpochFeat2_diff)
        diff_sum(i,j)=abs(EpochFeat2_diff(i,1,j))/EpochFeat1_mean(i,1)+...
            abs(EpochFeat2_diff(i,5,j))/EpochFeat1_mean(i,5)+...
            abs(EpochFeat2_diff(i,7,j))/EpochFeat1_mean(i,7)+...
            abs(EpochFeat2_diff(i,12,j))/EpochFeat1_mean(i,12);
        diff_sum(i,j)=diff_sum(i,j)/4;
        
        
        idx=[1+StEpochTime(i)*fs:(StEpochTime(i)+opt1.twin)*fs];
        
        
        
        
        
        %% apnea happens but no apean after the epoch. It should be deleted 
        if labelEpoch(i)==nan
            flag(i,j)==0;
        end
        
        
        
        %% no apnea epoch
        if labelEpoch(i)==0
            labelMark(1)=labelMark(1)+1;
            
            %              if
            %                     flag(i,j)=0;
            %                     NoMeetSelf(1)= NoMeetSelf(1)+1;
            
            
            
            if abs(EpochFeat2_diff(i,12,j))<opt1.psdDiff  &  ...
                    abs(EpochFeat2_diff(i,1,j))<opt1.p(1) & abs(EpochFeat2_diff(i,2,j))<opt1.p(2) &...
                    abs(EpochFeat2_diff(i,3,j))<opt1.p(1) & abs(EpochFeat2_diff(i,4,j))<opt1.p(2) &...
                    abs(EpochFeat2_diff(i,5,j))<opt1.p(5) & abs(EpochFeat2_diff(i,6,j))<opt1.p(6)&...
                    abs(EpochFeat2_diff(i,7,j))<opt1.p(5) & abs(EpochFeat2_diff(i,8,j))<opt1.p(6)&...
                    abs(EpochFeat2_diff(i,13,j))<opt1.p(7)
                
                
                
                flag(i,j)=1;
                MeetLimit(1)=MeetLimit(1)+1;
                
                if  EpochFeat2(i,15,j)>opt1.pself(5)  % covPP self thrshold
                    flag(i,j)=0;
                    NoMeetP(4)= NoMeetP(4)+1;
                    
                end
                
                
                if  EpochFeat2(i,12,j)<opt1.pself(4)  % psd self thrshold
                    flag(i,j)=0;
                    NoMeetP(3)= NoMeetP(3)+1;
                    
                end
                
                if   EpochFeat2(i,4,j)> opt1.pself(3)
                    flag(i,j)=0;
                    NoMeetP(2)= NoMeetP(2)+1;
                    
                end
                
                if   EpochFeat2(i,2,j)> opt1.pself(2)
                    flag(i,j)=0;
                    NoMeetP(1)= NoMeetP(1)+1;
                    
                end
                
                if  diff_sum(i,j)> opt1.diffTh
                    NoMeetDiff(2)= NoMeetDiff(1)+1;
                    flag(i,j)=0;
                end
              
                
                
            
            end
            
            % record number of t/f in no apnea epoches
            if flag(i,j)==0;
                falseNum(1)=falseNum(1)+1;
                
            else
                trueNum(1)=trueNum(1)+1;
                NCSch(j)=NCSch(j)+1;
            end
        end
        % apnea epoch
        
        if labelEpoch(i)== 1 | labelEpoch(i)== 2 | labelEpoch(i)== 3 
            labelMark(2)=labelMark(2)+1;
            if abs(EpochFeat2_diff(i,12,j))<opt2.psdDiff  &  ...
                    abs(EpochFeat2_diff(i,1,j))<opt2.p(1) & abs(EpochFeat2_diff(i,2,j))<opt2.p(2) &...
                    abs(EpochFeat2_diff(i,3,j))<opt2.p(1) & abs(EpochFeat2_diff(i,4,j))<opt2.p(2) &...
                    abs(EpochFeat2_diff(i,5,j))<opt2.p(5) & abs(EpochFeat2_diff(i,6,j))<opt2.p(6)&...
                    abs(EpochFeat2_diff(i,7,j))<opt2.p(5) & abs(EpochFeat2_diff(i,8,j))<opt2.p(6)&...
                    abs(EpochFeat2_diff(i,13,j))<opt2.p(7)
                
                
                flag(i,j)=1;
                MeetLimit(2)=MeetLimit(2)+1;
                if  diff_sum(i,j)> opt2.diffTh
                    NoMeetDiff(2)= NoMeetDiff(2)+1;
                    flag(i,j)=0;
                end
               
                if  EpochFeat2(i,12,j)<opt2.pself
                    flag(i,j)=0;
                    NoMeetP(2)= NoMeetP(2)+1;
                    
                end
                
                
                
                
            end
            % record number of t/f in apnea epoches
            if flag(i,j)==0;
                falseNum(2)=falseNum(2)+1;
                
            else
                trueNum(2)=trueNum(2)+1;
                NCSch(j)=NCSch(j)+1;
            end
            
        end
        %% 
          if labelEpoch(i)== 4 | labelEpoch(i)== 5 | labelEpoch(i)== 6
              labelMark(3)=labelMark(3)+1;
            if abs(EpochFeat2_diff(i,12,j))<opt3.psdDiff  &  ...
                    abs(EpochFeat2_diff(i,1,j))<opt3.p(1) & abs(EpochFeat2_diff(i,2,j))<opt3.p(2) &...
                    abs(EpochFeat2_diff(i,3,j))<opt3.p(1) & abs(EpochFeat2_diff(i,4,j))<opt3.p(2) &...
                    abs(EpochFeat2_diff(i,5,j))<opt3.p(5) & abs(EpochFeat2_diff(i,6,j))<opt3.p(6)&...
                    abs(EpochFeat2_diff(i,7,j))<opt3.p(5) & abs(EpochFeat2_diff(i,8,j))<opt3.p(6)&...
                    abs(EpochFeat2_diff(i,13,j))<opt3.p(7)
                
                
                flag(i,j)=1;
                MeetLimit(3)=MeetLimit(3)+1;
                if  diff_sum(i,j)> opt3.diffTh
                    NoMeetDiff(3)= NoMeetDiff(3)+1;
                    flag(i,j)=0;
                end
              
                
                %opt23.pself is only used for psd self threshold 
                  if  EpochFeat2(i,12,j)<opt3.pself
                    flag(i,j)=0;
                    NoMeetP(3)= NoMeetP(3)+1;
                    
                end
                
                
         
            end
            % record number of t/f in apnea epoches
            if flag(i,j)==0;
                falseNum(3)=falseNum(3)+1;
                
            else
                trueNum(3)=trueNum(3)+1;
                NCSch(j)=NCSch(j)+1;
            end
            
        end
        
    end
end

property.trueRate=length(flag(flag~=0))/size(flag,1)/size(flag,2);
property.trueRateNoApneaEpoch=trueNum(1)/(trueNum(1)+falseNum(1));
property.trueRateApneaEpoch(1)=trueNum(2)/(trueNum(2)+falseNum(2));  %apnea but not OSA CSA 
property.trueRateApneaEpoch(2)=trueNum(3)/(trueNum(3)+falseNum(3));  % apnea of OSA CSA 
property.NoMeetP1=NoMeetP;
property.MeetLimit=MeetLimit;
property.diff=diff_sum;
property.NoMeetDiff=NoMeetDiff;
property.labelMark=labelMark;
property.trueNCSch=NCSch./length(flag(flag==1));
[r,c]=find(flag==1);
property.trueNCS_UniqEpoch=length(unique(r));





end