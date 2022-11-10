function  [h]=PlotTrueTimeStamp(flag,labelEpoch,StEpochTime,validEpochIdx,opt)

idx0=find(validEpochIdx==0 );



    flag(idx0,:)=0;
    [r,c]=find(flag==1);
    t=1:StEpochTime(end)+opt.twin;
    tStamp=zeros(length(t),1);
    tlabel=zeros(length(t),2);
    for i=1:length(r)   % total number of true NCS epoch
     
        label(i)=labelEpoch(r(i));
        
        StTimeEpoch(i)=StEpochTime(r(i));
        
        
        tStamp(StTimeEpoch(i)+1:StTimeEpoch(i)+opt.twin)=1;
        if label(i)==0
            tlabel(StTimeEpoch(i)+1:StTimeEpoch(i)+opt.twin,1)=1;
        
        else 
            tlabel(StTimeEpoch(i)+1:StTimeEpoch(i)+opt.twin,2)=1;
            
        end
    
    end
  sz=10;  
  h(1)=figure();
  
  
  fill(t,tStamp,'r','EdgeColor','none');
  xlim([1 t(end)]);
  xlabel('time (s)','FontSize',sz);
  ylabel('True NCS Epoch','FontSize',sz);
  set(gcf,'Position',[100,101,1200,150]);
  
  
  r=length(tStamp(tStamp==1))/length(tStamp);
  txt1=['True NCS time rate:',num2str(r,2)];
  title(txt1,'FontSize',sz);
  
  h(2)=figure();
  

  fill(t,tlabel(:,1),'r','EdgeColor','none');
  hold on
  fill(t,tlabel(:,2),'g','EdgeColor','none');
  
  xlim([1 t(end)]);
  xlabel('time (s)','FontSize',sz);
  ylabel('True NCS Epoch','FontSize',sz);
  set(gcf,'Position',[100,100,1200,150]);
  legend('normal','disorder');
  
  r=length(tStamp(tStamp==1))/length(tStamp);
  txt1=['True NCS time rate:',num2str(r,2)];
  title(txt1,'FontSize',sz);

end