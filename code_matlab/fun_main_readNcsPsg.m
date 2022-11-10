%% predict the apnea in advance. 1~2 epoches before apnea happens


clear 
    %% read NCS and PSG data synchronize to mutual time duration 
  for i=1:30
   
 
 fileName1=['NCS0',num2str(i,'%02d')];
 fileName2=[num2str(i,'%02d')];
 main_process(fileName1,fileName2)

  end    
    
  
 for i=1:30
   
 fileName1=['NCS0',num2str(i,'%02d')];
 fileName2=[num2str(i,'%02d'),'_2'];
 main_process(fileName1,fileName2)

  end 

function main_process(fileName1,fileName2) 

  
close all 
 

DataVer='Feat_pred_v4';
DataVer='Feat_detect_v5';

%save features 
SavePath=['C:\Users\zz587-admin\Documents\sleep center\CNN_algorithm\Output\',DataVer,'\'];
% input NCS and PSG waveform. 
DataPath=['C:\Users\zz587-admin\Documents\sleep center\CNN_algorithm\NCS_PSG_data\',fileName2,'.mat'];

PredOpt=0; 
%PredOpt=3;%detect=0  pred =2,3 
if PredOpt==0
    opt.twin=40;opt.twinMove=30;
end
if PredOpt~=0 
    opt.twin=40;opt.twinMove=15;
end
opts1.smfilt=1;  % whehter to add smooth filter
fs=25;
if ~exist(DataPath,'file') 
    fprintf('not found NCS case %s\n',fileName2);
    return
    
end

load(DataPath); 
fprintf('case: %s \n',fileName2);
%%old version  of generation data each time. time consuming

% input sleep center data 
dataPath1=['C:\Users\zz587-admin\Documents\sleep center\data_reference\NCS EDF records\'];
dataPath2=['C:\Users\zz587-admin\Documents\sleep center\matFile\scpair\'];

filePathName1 = [dataPath1,fileName1,'.edf'];
filePathName2 = [dataPath2,fileName2,'.mat'];




[NCSdata,StTimeNCS] = readNCSmat(filePathName2,fs);% The order is [amp_tx1rx1 amp_tx2rx2 ph_tx1rx1 ph_tx2rx2]
[PSGdata,PSGspO2,StTimePSG,scorePSG,scorePSG_all,annotations]=readPSG(filePathName1,fs);% PSGdata len=t(sec)*fs ; score col 1: time sample  col 2:end : label 0-12  (0/1)
EndTimeNCS=StTimeNCS +seconds(size(NCSdata,1)/fs);
EndTimePSG=StTimePSG +seconds(size(PSGdata,1)/fs);

if EndTimeNCS<StTimePSG | EndTimePSG<StTimeNCS
   fprintf('no mutual time \n')
   return
end
PSGdata_all=[PSGdata PSGspO2];
[NCSdata,PSGdata_all,score,score_all,StTimeSyn]=synchNCS_PSG(NCSdata,StTimeNCS,PSGdata_all,StTimePSG,scorePSG,scorePSG_all,fs);

save(DataPath,'NCSdata','PSGdata_all','PSGdata','PSGspO2','score','score_all',...
    'StTimeSyn','annotations');

%% pre-process data with filter,  transfer into epoch data and extract features, label epoches

NCSdata_nofilt=NCSdata;
PSGdataCmp_nofilt=PSGdata_all(:,1:3);

opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 2; opts1.fstLP = 2.5; 
%opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
NCSdata = filterLpHp(NCSdata,fs,opts1);
PSGdata_all(:,1:3) = filterLpHp(PSGdata_all(:,1:3),fs,opts1);

% filter by smooth filter? Include noise? 

if opts1.smfilt==1
    NCSdata=sgolayfilt(NCSdata,4,61);
    PSGdata_all(:,1:3)=sgolayfilt(PSGdata_all(:,1:3),4,61);
end

% feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
%          skew_mean kurt_mean entro...
%            per_power cycle fmax CaliDel];
opt.psdBand=[0.05,0.5];
opt.minInterceptDist = 0.3; opt.tWinPk=4;% paramters in sub function :peakDet3
opt.calib=1;opt.calibT=[1 opt.twin];opt.calibMinPkRatio=0.4; %parameter for peak det. calibrate out p-p< ratio*mean
opt.lowTh=0.92;%opt.featNum=16;
opt.featNum=29;
PSGdataCmp=PSGdata_all(:,1:3);
PSGdataSpO2=PSGdata_all(:,4);


[EpochFeatPSGCmp,StEpochTimePSG,EpochDataPSG]=windowFeat(PSGdataCmp,fs,opt);
[EpochFeatPSGSpO2,~,EpochDataPSGSpO2]=windowFeat_SpO2(PSGdataSpO2,fs,opt);

[EpochFeatNCS,StEpochTimeNCS,EpochDataNCS]=windowFeat(NCSdata,fs,opt);


%% extra high frequency band features 

opt.HFband=[1 2;2 5; 5 8;8 12.5]; opt.featNumHF=8; 
[EpochFeatNCS_HF]=windowFeat_HF(NCSdata_nofilt,fs,opt);
[EpochFeatPSGCmp_HF]=windowFeat_HF(PSGdataCmp_nofilt,fs,opt);

%concate features 
EpochFeatNCS=cat(2,EpochFeatNCS,EpochFeatNCS_HF);
EpochFeatPSGCmp=cat(2,EpochFeatPSGCmp,EpochFeatPSGCmp_HF);
%% label epochs 

% opt1.twin=90;opt1.classMethod=1;opt1.labelRaio=0.25;
opt1.twin=opt.twin;opt1.labelRaio=0.25;opt1.classMethod=2;opt1.ratioThSnore=0.5;
[labelEpochOld,ratio,validEpochIdx]=EpochScoreLabel(score_all,score,StEpochTimePSG,opt1,fs);  % ratio =[apneaRaio(in epoch)  apeanRatio(for entire time snore include)];
%valid Epoch idx =0 when more than half score=0 


%% predict label. In advance of  #  epoches 
% PredOpt=4: Apnea epoch excluded. Only use normal epoch to predict In advance of 6 epoches   
% PredOpt=0  no prediction 
labelEpoch=EpochPredict(labelEpochOld,PredOpt);


 
%% Compare features and waveforms , and find True NCS epoches   plot example epoches

opt21.psdDiff=30; opt21.diffTh=0.4;opt21.twin=opt.twin;%opt21.pself=[60,6,0.22,65,60];
% opt21.pself=[60,6,0.22,65,60]; old version   p= self [meanBR  stdBR stdpp psd covPP] 

opt21.p=[20,8,1.5,1,2,1.5,8]; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 


opt22.psdDiff=40; opt22.diffTh=0.5;       % for label = 1 2 3
opt22.p=[30,8,2,1,2,1.5,10]; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt22.pself=50;

opt23.psdDiff=40; opt23.diffTh=0.5;       % for label = 4 5 6
opt23.p=[30,8,2,1,2,1.5,10]; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt23.pself=50;


% if detection: 
if PredOpt==0   
   opt21.pself=[60,6,0.22,70,0.3]; 
   opt22.pself=70;opt23.pself=70;
   opt23.psdDiff=25; opt22.psdDiff=25;
   opt23.p=[25,8,2,1,2,1.5,10];opt22.p=[25,8,2,1,2,1.5,10];
   fprintf('PredOpt==0');%% for prediction. change back to 50 for detection
end

% if prediction : 
if PredOpt~=0   
   opt21.pself=[60,5,0.2,80,0.3]; 

end

if PredOpt==5 % > 6 epoch, results will be worse
   opt21.pself=[60,8,0.3,60,0.35]; 
   opt22.pself=70;opt23.pself=70;
   opt23.psdDiff=30; opt22.psdDiff=30;
   opt23.p=[30,8,2,1,2,1.5,10];opt22.p=[30,8,2,1,2,1.5,10];
end




% old version of comparison 
%[flagNCS,property,PSGoptIdx]=featComp_ApneaLabel(EpochFeatPSGCmp,EpochFeatNCS,opt21,opt22,opt23,PSGdataCmp,NCSdata,StEpochTimeNCS,fs,labelEpoch);
%property: TrueRate: True NCS label rate  diff:mean of diff
[flagNCS,property,PSGoptIdx]=featComp_self(EpochFeatPSGCmp,EpochFeatNCS,opt21,opt22,opt23,StEpochTimeNCS,fs,labelEpoch);


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


%% Output NCS PSG labels and features  


opt4.ncs=1; 
[featureNCS_all,labelNCS_all]=outputResult(flagNCS,EpochFeatNCS,EpochFeatPSGSpO2,labelEpoch,opt4,validEpochIdx); 
% plot true NCS Epoch time stamps 
h0=PlotTrueTimeStamp(flagNCS,labelEpoch,StEpochTimeNCS,validEpochIdx,opt);


opt4.ncs=2;opt4.psgIdx=PSGoptIdx;
[featurePSG2_all,labelPSG2_all]=outputResult(flagNCS,EpochFeatPSGCmp,EpochFeatPSGSpO2,labelEpoch,opt4,validEpochIdx);

%feature PSG 2, the same epoches as collected from NCS
% output data: feature: epochNum *feature Num     label: epochNum (0,1)
opt4.ncs=0;opt4.psgIdx=PSGoptIdx;
[featurePSG_all,labelPSG_all]=outputResult(0,EpochFeatPSGCmp,EpochFeatPSGSpO2,labelEpoch,opt4,validEpochIdx);
TrueNCSepochApeanRatio=length(find(labelNCS_all~=0))/length(labelNCS_all);

fprintf('duration(h): %.2f, apneaRatio %.2f, true NCS epoch rate %.2f , TrueNCSepochApeanRatio %.2f \n',...
hours(seconds(size(PSGdata_all,1)/fs)),ratio.apneaRaioEpoch,property.trueRate,TrueNCSepochApeanRatio );



%% Output NCS PSG  waveform data directly   
%output 2 channels
%  ch1: waveform (norm to [-1,1])     
%  ch2: SpO2 *0.01 

opt4.ncs=1;
[NCS_data,labelNCS]=outputWaveform(flagNCS,EpochDataNCS,EpochDataPSGSpO2,labelEpoch,opt4,validEpochIdx); 

%feature PSG 2, the same epoches as collected from NCS
opt4.ncs=2;opt4.psgIdx=PSGoptIdx;
[PSG2_data,labelPSG2]=outputWaveform(flagNCS,EpochDataPSG,EpochDataPSGSpO2,labelEpoch,opt4,validEpochIdx); 

opt4.ncs=0;opt4.psgIdx=PSGoptIdx;
[PSG_data,labelPSG]=outputWaveform(flagNCS,EpochDataPSG,EpochDataPSGSpO2,labelEpoch,opt4,validEpochIdx); 



%% Sve data and figs


FeatureName=["mean(br)", "std(br)", "mean(pp)" ,"std(pp)" ,"mean(in)" ,"std(in)" ,"mean(ex)" ,'std(ex)',...
                'skew_mean' , 'kurt_mean' ,'entro',...
                'per_power', 'cycle', 'covBR' ,'covPP', 'void_t',...
                    'max(in)','max(ex)','max(br)','min(br)','min(pp)',...
                    'Cor(br)','SD(br)','Cor(pp)','SD(pp)','Cor(in)','SD(in)','Cor(ex)','SD(ex)',...
                    'psd1','psd2','psd3','psd4','fmax1','fmax2','fmax3','fmax4',...
                    'mean(spo2)' ,'std(spo2)', 'percent(spo2)', 'min(spo2)'];

SaveProp.FeatureName=FeatureName;
SaveProp.SartTime=StTimeSyn;SaveProp.Duration=hours(seconds(size(PSGdata,1)/fs));
SaveProp.ApneaRaio=ratio; SaveProp.featCompare=property;SaveProp.apneaLabelRaio=opt1.labelRaio;
SaveProp.apneaLabelRaio=opt1.labelRaio;SaveProp.classMethod=opt1.classMethod;
SaveProp.TrueNCSepochOpt={opt21,opt22,opt23};SaveProp.TrueNCSepochApeanRatio=TrueNCSepochApeanRatio;
SaveProp.trueNCSRate=property.trueRate;SaveProp.flagNCS=flagNCS;SaveProp.Filter=opts1;

%include all feature extracted output = case num * feature num (20)
Feature.featureNCS=featureNCS_all;Feature.labelNCS=labelNCS_all;
Feature.featurePSG2=featurePSG2_all;Feature.labelPSG2=labelPSG2_all;
Feature.featurePSG=featurePSG_all;Feature.labelPSG=labelPSG_all;



save([SavePath,'segments\',fileName2,'.mat'],'SaveProp','Feature',...
      'NCS_data','labelNCS',...
      'PSG2_data','labelPSG2','PSG_data','labelPSG');

figName1 = [SavePath,'fig\',fileName2,'waveform_0'];
print(h1,[figName1,'.tiff'],'-dtiff','-r300');
savefig(h1,[figName1,'.fig']);
if h2==0
else
figName2 = [SavePath,'fig\',fileName2,'waveform_1'];
print(h2,[figName2,'.tiff'],'-dtiff','-r300');
savefig(h2,[figName2,'.fig']);
end
if h3==0
else
figName3 = [SavePath,'fig\',fileName2,'waveform_2'];
print(h3,[figName3,'.tiff'],'-dtiff','-r300');
savefig(h3,[figName3,'.fig']);
end
if h4==0
else
figName4 = [SavePath,'fig\',fileName2,'waveform_3'];
print(h4,[figName4,'.tiff'],'-dtiff','-r300');
savefig(h4,[figName4,'.fig']);
end

figName01 = [SavePath,'fig\timeStamp\',fileName2,'timeStamp_1'];
print(h0(1),[figName01,'.tiff'],'-dtiff','-r300');
savefig(h0(1),[figName01,'.fig']);
figName02 = [SavePath,'fig\timeStamp\',fileName2,'timeStamp_2'];
print(h0(2),[figName02,'.tiff'],'-dtiff','-r300');
savefig(h0(2),[figName02,'.fig']);

%clear
  end