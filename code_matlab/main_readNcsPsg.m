% for i=1:30
% 
%    
% fileName1=['NCS0',num2str(i,'%02d')];
% fileName2=num2str(i,'%02d');
% 
% main_fun(fileName1,fileName2);
% end
% 
% 
% function main_fun(fileName1,fileName2)

    %% read NCS and PSG data synchronize to mutual time duration 
   

 for i=2 
 fileName1=['NCS0',num2str(i,'%02d')];
 fileName2=[num2str(i,'%02d')];
 
    
dataPath1=['C:\Users\zz587-admin\Documents\sleep center\data_reference\NCS EDF records\'];
dataPath2=['C:\Users\zz587-admin\Documents\sleep center\matFile\scpair\'];

filePathName1 = [dataPath1,fileName1,'.edf'];
filePathName2 = [dataPath2,fileName2,'.mat'];

if ~exist(filePathName2,'file') 
    fprintf('not found NCS case %s\n',fileName2);
    continue
    
end
fs=20;%PSG airflow has best quality, chest and abdomen contains noise
[NCSdata,StTimeNCS] = readNCSmat(filePathName2,fs);% The order is [amp_tx1rx1 amp_tx2rx2 ph_tx1rx1 ph_tx2rx2]
[PSGdata,StTimePSG,scorePSG,annotations]=readPSG(filePathName1,fs);% PSGdata len=t(sec)*fs ; score col1: time(sec) col2:label 0-10
EndTimeNCS=StTimeNCS +seconds(size(NCSdata,1)/fs);
EndTimePSG=StTimePSG +seconds(size(PSGdata,1)/fs);
fprintf('case: %s \n',fileName2);
if EndTimeNCS<StTimePSG | EndTimePSG<StTimeNCS
   fprintf('no mutual time \n')
   return
end
[NCSdata,PSGdata,score,StTimeSyn]=synchNCS_PSG(NCSdata,StTimeNCS,PSGdata,StTimePSG,scorePSG,fs);
%% pre-process data with filter,  transfer into epoch data and extract features, label epoches


opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
NCSdata = filterLpHp(NCSdata,fs,opts1);
PSGdata = filterLpHp(PSGdata,fs,opts1);
NCSdata_nofilt=NCSdata;
% filter by smooth filter? Include noise? 
NCSdata=sgolayfilt(NCSdata,4,61);
PSGdata=sgolayfilt(PSGdata,4,61);

% feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
%          skew_mean kurt_mean entro...
%            per_power cycle fmax CaliDel];
opt.twin=90;opt.twinMove=30;opt.psdBand=[0.05,0.5];
opt.minInterceptDist = 0.3; opt.tWinPk=4;% paramters in sub function :peakDet3
opt.calib=1;opt.calibT=[1 90];opt.calibMinPkRatio=0.4;  %parameter for peak det. calibrate out p-p< ratio*mean
[EpochFeatPSG,StEpochTimePSG,EpochDataPSG]=windowFeat(PSGdata,fs,opt);
[EpochFeatNCS,StEpochTimeNCS,EpochDataNCS]=windowFeat(NCSdata,fs,opt);

%[EpochFeatNCS_nofilt,~,~]=windowFeat(NCSdata_nofilt,fs,opt);


% opt1.twin=90;opt1.classMethod=1;opt1.labelRaio=0.25;
opt1.twin=90;opt1.labelRaio=0.2;opt1.classMethod=2;opt1.snoreWeight=0.2;
[labelEpoch,ratio,validEpochIdx]=EpochScoreLabel(score,StEpochTimePSG,opt1);  % ratio =[apneaRaio(in epoch)  apeanRatio(for entire time snore include)];
%valid Epoch idx =0 when more than half score=0 
%% Compare features and waveforms , and find True NCS epoches   plot example epoches

opt21.psdDiff=30; opt21.diffTh=0.4;opt21.twin=90;opt21.pself=[60,7,0.28,65,60];% p= self [meanBR  stdBR covpp psd] 
opt21.p=[15,6,1,0.12,1.2,1.2,6];opt21.corTh=0.13; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt21.dtwTh=400;

opt22.psdDiff=30; opt22.diffTh=0.4;
opt22.p=[15,5,2,0.3,2,2,10];opt22.corTh=0.08; % p=diff in [meanBR  stdBR meanpp stdpp meanIn stdIn  cycle] 
opt22.dtwTh=500;
%[flagNCS,property,PSGoptIdx]=featComp(EpochFeatPSG,EpochFeatNCS,opt21,PSGdata,NCSdata,StEpochTimeNCS,fs); 
[flagNCS,property,PSGoptIdx]=featComp_ApneaLabel(EpochFeatPSG,EpochFeatNCS,opt21,opt22,PSGdata,NCSdata,StEpochTimeNCS,fs,labelEpoch);
%property: TrueRate: True NCS label rate  diff:mean of diff



opt3=opt;
opt3.seqPlot=0.3;  
%seqPlot = plot index in epoches     set a threshold to diff.
 %parameter for peak det. calibrate out p-p< ratio*mean
opt3.diffTh=0.4; opt3.dtwTh=600; opt3.corTh=0.1; opt3.snoreWeight=opt1.snoreWeight;
h=plotEpoch(PSGdata,NCSdata,StEpochTimeNCS,score,opt3,fs,labelEpoch,flagNCS,property,PSGoptIdx,EpochFeatNCS);
% opt4=opt3;opt4.seqPlot=0.7;
% h1=plotEpoch(PSGdata,NCSdata,StEpochTimeNCS,score,opt4,fs,labelEpoch,flagNCS,property,PSGoptIdx);
%% Output NCS PSG labels and features  save results

opt4.ncs=1;
[featureNCS_all,labelNCS_all]=outputResult(flagNCS,EpochFeatNCS,labelEpoch,opt4,validEpochIdx); 
% output data: feature: epochNum *feature Num     label: epochNum (0,1)
opt4.ncs=0;opt4.psgIdx=PSGoptIdx;
[featurePSG_all,labelPSG_all]=outputResult(0,EpochFeatPSG,labelEpoch,opt4,validEpochIdx);
TrueNCSepochApeanRatio=length(find(labelNCS_all==1))/length(labelNCS_all);

fprintf('duration(h): %.2f, apneaRatio %.2f, true NCS epoch rate %.2f , TrueNCSepochApeanRatio %.2f \n',...
hours(seconds(size(PSGdata,1)/fs)),ratio.apneaRaioEpoch,property.trueRate,TrueNCSepochApeanRatio );


SaveProp.SartTime=StTimeSyn;SaveProp.Duration=hours(seconds(size(PSGdata,1)/fs));
SaveProp.ApneaRaio=ratio; SaveProp.featCompare=property;SaveProp.apneaLabelRaio=opt1.labelRaio;
SaveProp.apneaLabelRaio=opt1.labelRaio;SaveProp.classMethod=opt1.classMethod;
SaveProp.TrueNCSepochOpt=opt21;SaveProp.TrueNCSepochApeanRatio=TrueNCSepochApeanRatio;
SaveProp.trueNCSRate=property.trueRate;SaveProp.flagNCS=flagNCS;


SavePath='C:\Users\zz587-admin\Documents\sleep center\OutputFeat\segments\';
save([SavePath,fileName2,'.mat'],'SaveProp','featureNCS_all','labelNCS_all','featurePSG_all','labelPSG_all');
figName1 = [SavePath,'fig\',fileName2,'waveform'];
print(h(1),[figName1,'.tiff'],'-dtiff','-r300');
savefig(h(1),[figName1,'.fig']);

figName2 = [SavePath,'fig\',fileName2,'peakDet'];
print(h(2),[figName2,'.tiff'],'-dtiff','-r300');
savefig(h(2),[figName2,'.fig']);

clear
 end