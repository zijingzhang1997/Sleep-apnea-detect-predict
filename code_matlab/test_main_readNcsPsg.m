%% predict the apnea in advance. 1~2 epoches before apnea happens
%opt.twin=30;opt.twinMove=10;
    %% read NCS and PSG data synchronize to mutual time duration 
 
    
    

 for i=1:30
     close all 
 fileName1=['NCS0',num2str(i,'%02d')];
 fileName2=[num2str(i,'%02d')];
fileName2=[num2str(i,'%02d'),'_2'];
    
dataPath1=['C:\Users\zz587-admin\Documents\sleep center\data_reference\NCS EDF records\'];
dataPath2=['C:\Users\zz587-admin\Documents\sleep center\matFile\scpair\'];

filePathName1 = [dataPath1,fileName1,'.edf'];
filePathName2 = [dataPath2,fileName2,'.mat'];

if ~exist(filePathName2,'file') 
    fprintf('not found NCS case %s\n',fileName2);
    continue
    
end


fs=25;%PSG airflow has best quality, chest and abdomen contains noise
[NCSdata,StTimeNCS] = readNCSmat(filePathName2,fs);% The order is [amp_tx1rx1 amp_tx2rx2 ph_tx1rx1 ph_tx2rx2]
[PSGdata,PSGspO2,StTimePSG,scorePSG,scorePSG_all,annotations]=readPSG(filePathName1,fs);% PSGdata len=t(sec)*fs ; score col 1: time sample  col 2:end : label 0-12  (0/1)
EndTimeNCS=StTimeNCS +seconds(size(NCSdata,1)/fs);
EndTimePSG=StTimePSG +seconds(size(PSGdata,1)/fs);
fprintf('case: %s \n',fileName2);
if EndTimeNCS<StTimePSG | EndTimePSG<StTimeNCS
   fprintf('no mutual time \n')
   return
end
PSGdata_all=[PSGdata PSGspO2];
[NCSdata,PSGdata_all,score,score_all,StTimeSyn]=synchNCS_PSG(NCSdata,StTimeNCS,PSGdata_all,StTimePSG,scorePSG,scorePSG_all,fs);
%% pre-process data with filter,  transfer into epoch data and extract features, label epoches


opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 2; opts1.fstLP = 2.5;
NCSdata = filterLpHp(NCSdata,fs,opts1);
PSGdata_all(:,1:3) = filterLpHp(PSGdata_all(:,1:3),fs,opts1);
NCSdata_nofilt=NCSdata;
% filter by smooth filter? Include noise? 
NCSdata=sgolayfilt(NCSdata,4,61);
PSGdata_all(:,1:3)=sgolayfilt(PSGdata_all(:,1:3),4,61);

% feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
%          skew_mean kurt_mean entro...
%            per_power cycle fmax CaliDel];
opt.twin=40;opt.twinMove=15;opt.psdBand=[0.05,0.5];
opt.minInterceptDist = 0.3; opt.tWinPk=4;% paramters in sub function :peakDet3
opt.calib=1;opt.calibT=[1 opt.twin];opt.calibMinPkRatio=0.4; %parameter for peak det. calibrate out p-p< ratio*mean
opt.lowTh=0.92;opt.featNum=16;
PSGdataCmp=PSGdata_all(:,1:3);
PSGdataSpO2=PSGdata_all(:,4);
[EpochFeatPSGCmp,StEpochTimePSG,EpochDataPSG]=windowFeat(PSGdataCmp,fs,opt);
[EpochFeatPSGSpO2,~,~]=windowFeat_SpO2(PSGdataSpO2,fs,opt);



[EpochFeatNCS,StEpochTimeNCS,EpochDataNCS]=windowFeat(NCSdata,fs,opt);

%[EpochFeatNCS_nofilt,~,~]=windowFeat(NCSdata_nofilt,fs,opt);


% opt1.twin=90;opt1.classMethod=1;opt1.labelRaio=0.25;
opt1.twin=opt.twin;opt1.labelRaio=0.25;opt1.classMethod=2;opt1.ratioThSnore=0.5;
[labelEpochOld,ratio,validEpochIdx]=EpochScoreLabel(score_all,score,StEpochTimePSG,opt1,fs);  % ratio =[apneaRaio(in epoch)  apeanRatio(for entire time snore include)];
%valid Epoch idx =0 when more than half score=0 


%% predict label. In advance of 2 epoches 
PredOpt=2;
labelEpoch=EpochPredict(labelEpochOld,PredOpt);



%% Compare features and waveforms , and find True NCS epoches   plot example epoches

opt21.psdDiff=30; opt21.diffTh=0.4;opt21.twin=opt.twin;opt21.pself=[60,6,0.22,75,60];
%old opt21.pself=[60,6,0.22,65,60]; p= self [meanBR  stdBR covpp psd_per] 
opt21.p=[15,6,1,0.12,1.2,1.2,4];opt21.corTh=0.13; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt21.dtwTh=400;

opt22.psdDiff=30; opt22.diffTh=0.4;       % for label = 1 2 3
opt22.p=[15,5,1,0.5,2,2,4];opt22.corTh=0; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt22.dtwTh=500;opt22.pself=50;

opt23.psdDiff=40; opt23.diffTh=0.5;       % for label = 4 5 6
opt23.p=[30,10,2,1,2,2,10];opt23.corTh=0; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt23.dtwTh=700;opt23.pself=50;
%[flagNCS,property,PSGoptIdx]=featComp(EpochFeatPSG,EpochFeatNCS,opt21,PSGdata,NCSdata,StEpochTimeNCS,fs); 
[flagNCS,property,PSGoptIdx]=featComp_ApneaLabel(EpochFeatPSGCmp,EpochFeatNCS,opt21,opt22,opt23,PSGdataCmp,NCSdata,StEpochTimeNCS,fs,labelEpoch);
%property: TrueRate: True NCS label rate  diff:mean of diff



opt3=opt;opt3.seqPlot=0.5;  
opt3.selectLabel=0;
%seqPlot = plot index in epoches   selectLabel
h1=plotEpoch(PSGdataCmp,NCSdata,PSGdataSpO2,StEpochTimeNCS,score,opt3,fs,labelEpoch,flagNCS,property,PSGoptIdx,EpochFeatNCS,fileName2);
opt3.selectLabel=1;
h2=plotEpoch(PSGdataCmp,NCSdata,PSGdataSpO2,StEpochTimeNCS,score,opt3,fs,labelEpoch,flagNCS,property,PSGoptIdx,EpochFeatNCS,fileName2);
opt3.selectLabel=2;
h3=plotEpoch(PSGdataCmp,NCSdata,PSGdataSpO2,StEpochTimeNCS,score,opt3,fs,labelEpoch,flagNCS,property,PSGoptIdx,EpochFeatNCS,fileName2);
opt3.selectLabel=3;
h4=plotEpoch(PSGdataCmp,NCSdata,PSGdataSpO2,StEpochTimeNCS,score,opt3,fs,labelEpoch,flagNCS,property,PSGoptIdx,EpochFeatNCS,fileName2);
%% Output NCS PSG labels and features  save results

opt4.ncs=1;
[featureNCS_all,labelNCS_all]=outputResult(flagNCS,EpochFeatNCS,EpochFeatPSGSpO2,labelEpoch,opt4,validEpochIdx); 
opt4.ncs=2;opt4.psgIdx=PSGoptIdx;
[featurePSG2_all,labelPSG2_all]=outputResult(flagNCS,EpochFeatPSGCmp,EpochFeatPSGSpO2,labelEpoch,opt4,validEpochIdx);
%feature PSG 2, the same epoches as collected from NCS
% output data: feature: epochNum *feature Num     label: epochNum (0,1)
opt4.ncs=0;opt4.psgIdx=PSGoptIdx;
[featurePSG_all,labelPSG_all]=outputResult(0,EpochFeatPSGCmp,EpochFeatPSGSpO2,labelEpoch,opt4,validEpochIdx);
TrueNCSepochApeanRatio=length(find(labelNCS_all~=0))/length(labelNCS_all);

fprintf('duration(h): %.2f, apneaRatio %.2f, true NCS epoch rate %.2f , TrueNCSepochApeanRatio %.2f \n',...
hours(seconds(size(PSGdata_all,1)/fs)),ratio.apneaRaioEpoch,property.trueRate,TrueNCSepochApeanRatio );


SaveProp.SartTime=StTimeSyn;SaveProp.Duration=hours(seconds(size(PSGdata,1)/fs));
SaveProp.ApneaRaio=ratio; SaveProp.featCompare=property;SaveProp.apneaLabelRaio=opt1.labelRaio;
SaveProp.apneaLabelRaio=opt1.labelRaio;SaveProp.classMethod=opt1.classMethod;
SaveProp.TrueNCSepochOpt=opt21;SaveProp.TrueNCSepochApeanRatio=TrueNCSepochApeanRatio;
SaveProp.trueNCSRate=property.trueRate;SaveProp.flagNCS=flagNCS;


SavePath='C:\Users\zz587-admin\Documents\sleep center\result\OutputFeat_predv2\';
save([SavePath,'segements\',fileName2,'.mat'],'SaveProp','featureNCS_all','labelNCS_all',...
      'featurePSG_all','labelPSG_all','featurePSG2_all','labelPSG2_all');
figName1 = [SavePath,'fig\',fileName2,'waveform_0'];
print(h1,[figName1,'.tiff'],'-dtiff','-r300');
%savefig(h(1),[figName1,'.fig']);
if h2==0
else
figName2 = [SavePath,'fig\',fileName2,'waveform_1'];
print(h2,[figName2,'.tiff'],'-dtiff','-r300');
end
if h3==0
else
figName3 = [SavePath,'fig\',fileName2,'waveform_2'];
print(h3,[figName3,'.tiff'],'-dtiff','-r300');
end
if h4==0
else
figName4 = [SavePath,'fig\',fileName2,'waveform_3'];
print(h4,[figName4,'.tiff'],'-dtiff','-r300');
end


clear
 end