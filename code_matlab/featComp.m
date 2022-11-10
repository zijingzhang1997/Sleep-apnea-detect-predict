function [flag,property,optIdx] =featComp(EpochFeat1,EpochFeat2,opt,PSGdata,NCSdata,StEpochTime,fs)
%feat 1 reference   feat 2 compare object
%feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
%          skew_mean kurt_mean entro...
%          per_power fmax];  compare [mean(br)1 mean(in)5 mean(ex)7  entro per_power12 ]

% EpochFeat1_mean=mean(EpochFeat1,3);

flag=zeros(size(EpochFeat2,1,3));
MeetLimit=0;NoMeetDiff=0;NoMeetCor=0;NoMeetDist=0;
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
           
         
     idx=[1+StEpochTime(i)*fs:(StEpochTime(i)+opt.twin)*fs];
     
     EpochTempNCS=NCSdata(idx,j);
     EpochTempPSG=PSGdata(idx,optIdx(i));
     %rescale to [0,1] in every epoch
     EpochTempNCS = rescale(EpochTempNCS,-1,1);
     EpochTempPSG = rescale(EpochTempPSG,-1,1);
     [c,~] = xcorr(EpochTempPSG,EpochTempNCS,fs*10,'normalized');
     cor(i,j)=max(c);
         
         
      if isfield(opt,'dtwTh')
          dist(i,j) = dtw(EpochTempNCS,EpochTempPSG); 
          
          
          
      end
      if abs(EpochFeat2_diff(i,12,j))<opt.psdDiff  &  abs(EpochFeat2_diff(i,1,j))<opt.p(1) & ...
           abs(EpochFeat2_diff(i,2,j))<opt.p(2) & abs(EpochFeat2_diff(i,5,j))<opt.p(3)&...
           abs(EpochFeat2_diff(i,6,j))<opt.p(4) & abs(EpochFeat2_diff(i,7,j))<opt.p(3)&...
           abs(EpochFeat2_diff(i,8,j))<opt.p(4)& abs(EpochFeat2_diff(i,13,j))<opt.p(5)
           flag(i,j)=1;
         MeetLimit=MeetLimit+1;
           if  diff_sum(i,j)> opt.diffTh
           NoMeetDiff= NoMeetDiff+1;
               flag(i,j)=0;
            end
      
         
         
   
             
         
         if cor(i,j)< opt.corTh
             NoMeetCor= NoMeetCor+1;
             flag(i,j)=0;
         end
         
           if isfield(opt,'dtwTh')
            if dist(i,j)>opt.dtwTh
              NoMeetDist= NoMeetDist+1;
             flag(i,j)=0;
            end
          end
         
         
         
         
       end   
         
  
  end
end

property.trueRate=length(flag(flag==1))/size(flag,1)/size(flag,2);

property.MeetLimit=MeetLimit;
property.diff=diff_sum;
property.NoMeetDiff=NoMeetDiff;
if isfield(opt,'corTh')
property.cor=cor;
property.NoMeetCor=NoMeetCor;
end
if isfield(opt,'dtwTh')
property.dist=dist;
property.NoMeetDist=NoMeetDist;
end



end