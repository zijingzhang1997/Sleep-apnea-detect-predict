function h=plotEpoch(PSGdata,NCSdata,PSGdata2,StEpochTime,score,opt,fs,labelEpoch,flagNCS,property,PSGoptIdx,EpochFeatNCS,fileName)
 
twin=opt.twin;
diff=property.diff;
num=1;
r=[];
c=[];
predTime=2;  % predict in advance of two epoches 
for i=1:size(flagNCS,1)
    for j=1:size(flagNCS,2)
        
        if opt.selectLabel==0
         if flagNCS(i,j)==1  & labelEpoch(i)==0
             r(num)=i;
             c(num)=j;
             num=num+1;
            
         end
        end
        if opt.selectLabel==1
         if flagNCS(i,j)==1  & (labelEpoch(i)==1 | labelEpoch(i)==2 | labelEpoch(i)==3)
             r(num)=i;
             c(num)=j;
             num=num+1;
            
         end 
        end
        if opt.selectLabel==2
         if flagNCS(i,j)==1  & (labelEpoch(i)==4 | labelEpoch(i)==5 )
             r(num)=i;
             c(num)=j;
             num=num+1;
             
         end 
        end
        if opt.selectLabel==3
         if flagNCS(i,j)==1  & ( labelEpoch(i)==6)
             r(num)=i;
             c(num)=j;
             num=num+1;
             
         end 
        end
    end
end

if num==1
    h=0;
    fprintf('no epoch found in %d label',opt.selectLabel);
    return
end
     
i=round(opt.seqPlot*length(r));
 %plot epoch result  identified by index 
idx=[1+StEpochTime(r(i))*fs:(StEpochTime(r(i))+twin)*fs];    %not sure about index



idx_pred=[1+StEpochTime(r(i))*fs:(StEpochTime(r(i))+twin*(predTime+1))*fs];
scoreEpoch=score(idx_pred,2);
scoreEpoch1=scoreEpoch(1:twin*fs);
scoreEpoch2=scoreEpoch(twin*fs+1:twin*fs*2);
scoreEpoch3=scoreEpoch(twin*fs*2+1:twin*fs*3);
apneaRate1=length(find(scoreEpoch1>=5))/length(scoreEpoch1);
apneaRate2=length(find(scoreEpoch2>=5))/length(scoreEpoch2);
apneaRate3=length(find(scoreEpoch3>=5))/length(scoreEpoch3);





SrTimeEpoch=hours(seconds(StEpochTime(r(i))));

LabelNamestr=["no apnea","snore","RERA+arousal","hypopnea","OSA","MA","CSA"];
LabelName=LabelNamestr((labelEpoch(r(i)))+1);

EpochfeatNCStemp=EpochFeatNCS(r(i),:,c(i));


EpochTempNCS=NCSdata(idx,c(i));
EpochTempPSG=PSGdata(idx,:);

EpochTempPSG2=PSGdata2(idx,:);
EpochTempNCS = rescale(EpochTempNCS,-1,1); 

t=((0:(length(EpochTempNCS)-1))/fs)';


EpochTempPSGopt=EpochTempPSG(:,PSGoptIdx(r(i)));
EpochTempPSGopt = rescale(EpochTempPSGopt,-1,1);
[featPSG,pkidxPSG]=PlotMaxMin_FeatInEpoch(EpochTempPSGopt,fs,opt);
[featNCS,pkidxNCS]=PlotMaxMin_FeatInEpoch(EpochTempNCS,fs,opt);

ncsChName=["amp Th","amp Ab","ph Th","ph Th"];
psgChName=["Airflow","chest","abdomen"];
h(1)=figure;
nfig=7;
sz=9;


subplot(nfig,1,3);
plot(t,EpochTempPSG(:,1),'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Airflow PSG','FontSize',sz)

subplot(nfig,1,4);
plot(t,EpochTempPSG(:,2),'LineWidth',0.5,'color','m');
xlabel('time (s)','FontSize',sz)
ylabel('Thorax PSG ','FontSize',sz)


subplot(nfig,1,5);
plot(t,EpochTempPSG(:,end),'LineWidth',0.5,'color','c');
xlabel('time (s)','FontSize',sz)
ylabel('Abdomen PSG','FontSize',sz)

subplot(nfig,1,6);
% spo2Baseline=EpochTempPSG2;
% spo2Baseline(:)=90;
plot(t,EpochTempPSG2,'LineWidth',0.5,'color','#D95319');
% hold on;
% plot(t,spo2Baseline,'LineStyle',':','LineWidth',1,'color','black');
xlabel('time (s)','FontSize',sz)
ylabel('SpO2 from PSG','FontSize',sz)


subplot(nfig,1,7);
t_pred=((0:(length(scoreEpoch)-1))/fs)';
scoreEpochBaseline=scoreEpoch;
scoreEpochBaseline(:)=4.5;
plot(t_pred,scoreEpoch,'LineWidth',1,'color','blue');
hold on
plot(t_pred,scoreEpochBaseline,'LineStyle',':','LineWidth',1,'color','black');
xlabel('time (s)','FontSize',sz)
ylabel('Score','FontSize',sz)
ylim([0 12])
title(['apnea rate Now: ',num2str(apneaRate1,3),'apnea rate Next: ',num2str(apneaRate2,3),'apnea rate Next 2: ',num2str(apneaRate3,3)],'FontSize',sz);
txt0=['case:',fileName,'epoch label: ',num2str(labelEpoch(r(i))),LabelName,' Start time(h): ',num2str(SrTimeEpoch,3)];




subplot(nfig,1,1);
plot(t, EpochTempNCS,'LineWidth',0.5,'color','red');
hold on;
pkMax=pkidxNCS.max;pkMin=pkidxNCS.min;
plot(t(pkMax),EpochTempNCS(pkMax),'^',...
     t(pkMin),EpochTempNCS(pkMin),'v');
xlabel('time (s)','FontSize',sz)
ylabel(['NCS '],'FontSize',sz)
txt1=['BR(BMP):',num2str(featNCS(1),3),', std BR',num2str(featNCS(2),3),', mean pp',num2str(featNCS(3),3),', std pp',num2str(featNCS(4),3),', psdRatio:',num2str(featNCS(12),4),', Ch:',ncsChName(c(i))];
title(join(txt1),'FontSize',sz);
set(gca,'fontsize', sz)

subplot(nfig,1,2);
plot(t,EpochTempPSGopt,'LineWidth',0.5,'color','green');
hold on;
pkMax=pkidxPSG.max;pkMin=pkidxPSG.min;
plot(t(pkMax),EpochTempPSGopt(pkMax),'^',...
     t(pkMin),EpochTempPSGopt(pkMin),'v');
xlabel('time (s)','FontSize',sz)
ylabel('PSG optimal','FontSize',sz)
txt2=['BR(BMP):',num2str(featPSG(1),3),', std BR:',num2str(featPSG(2),3),', mean pp:',num2str(featPSG(3),3),', std pp:',num2str(featPSG(4),3),', psdRatio:',num2str(featPSG(12),4),', Ch:',psgChName(PSGoptIdx(r(i)))];
title(join(txt2),'FontSize',sz)

txt3=['diff:', num2str(diff(r(i),c(i)),3)];

sgtitle([join(txt0) join(txt3)],'fontsize', sz);
set(gcf,'Position',[100,10,800,900]);
set(gca,'fontsize', sz);





end