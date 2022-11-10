Path='C:\Users\zz587-admin\Documents\sleep center\CNN_algorithm\Output\';
%FolderName='CNNFeat_pred_v4'; 

FolderName='Feat_detect_v4';% detection folder 

Path=[Path,FolderName,'\'];

dataPath=[Path,'segments\'];
SavePath=[Path,'featureAll\'];
SaveName=[FolderName,'_feature_data'];




filename=dir(dataPath);
labelNCS_all=[];
caseNumNCS=[];
labelPSG_all=[];
labelPSG2_all=[];
caseNumPSG=[];
caseNumPSG2=[];
featureNCS_all=[];
featurePSG_all=[];
featurePSG2_all=[];
NCSuniqEpoch=zeros(30,1);
num=0;
opt=1;
for i = 1:30
    fileName=[num2str(i,'%02d')];
    filePathName=[dataPath,fileName,'.mat'];
    
    if ~exist(filePathName,'file')
        fprintf('not found NCS case %s\n',fileName);
        continue
    end
    load(filePathName);
    
    featureNCS_seg=Feature.featureNCS;labelNCS_seg=Feature.labelNCS;
    featurePSG_seg=Feature.featurePSG;labelPSG_seg=Feature.labelPSG;
    featurePSG2_seg=Feature.featurePSG2;labelPSG2_seg=Feature.labelPSG2;
        
    NCSuniqEpoch(i)=SaveProp.featCompare.trueNCS_UniqEpoch;  
    
    if num==0
    labelNCS_all=labelNCS_seg;
    featureNCS_all=featureNCS_seg;
    labelPSG_all=labelPSG_seg;
    labelPSG2_all=labelPSG2_seg;
    featurePSG_all=featurePSG_seg;
    featurePSG2_all=featurePSG2_seg;
    caseNumNCS=zeros(length(labelNCS),1);
    caseNumNCS(:)=i;
    caseNumPSG=zeros(length(labelPSG),1);
    caseNumPSG(:)=i;
    caseNumPSG2=zeros(length(labelPSG2),1);
    caseNumPSG2(:)=i;
    else
        labelNCS_all=[labelNCS_all;labelNCS_seg];
        featureNCS_all=[featureNCS_all;featureNCS_seg];
        labelPSG_all=[labelPSG_all;labelPSG_seg];
        featurePSG_all=[featurePSG_all;featurePSG_seg];
        labelPSG2_all=[labelPSG2_all;labelPSG2_seg];
        featurePSG2_all=[featurePSG2_all;featurePSG2_seg];
        
        
        caseNumNCS=[caseNumNCS; ones(length(labelNCS_seg),1).*i];
        caseNumPSG=[caseNumPSG; ones(length(labelPSG_seg),1).*i];
        caseNumPSG2=[caseNumPSG2; ones(length(labelPSG2_seg),1).*i];
        
    end
    
    num=num+1;
    
end

for i = 1:30
    fileName=[num2str(i,'%02d'),'_2'];
    filePathName=[dataPath,fileName,'.mat'];
    
    if ~exist(filePathName,'file')
        
        continue
    end
    fprintf(' found NCS case %s\n',fileName);
    load(filePathName);
    featureNCS_seg=Feature.featureNCS;labelNCS_seg=Feature.labelNCS;
    featurePSG_seg=Feature.featurePSG;labelPSG_seg=Feature.labelPSG;
    featurePSG2_seg=Feature.featurePSG2;labelPSG2_seg=Feature.labelPSG2;
        
    NCSuniqEpoch(i)=SaveProp.featCompare.trueNCS_UniqEpoch;  
    
    if num==0
    labelNCS_all=labelNCS_seg;
    featureNCS_all=featureNCS_seg;
    labelPSG_all=labelPSG_seg;
    labelPSG2_all=labelPSG2_seg;
    featurePSG_all=featurePSG_seg;
    featurePSG2_all=featurePSG2_seg;
    caseNumNCS=zeros(length(labelNCS),1);
    caseNumNCS(:)=i;
    caseNumPSG=zeros(length(labelPSG),1);
    caseNumPSG(:)=i;
    caseNumPSG2=zeros(length(labelPSG2),1);
    caseNumPSG2(:)=i;
    else
        labelNCS_all=[labelNCS_all;labelNCS_seg];
        featureNCS_all=[featureNCS_all;featureNCS_seg];
        labelPSG_all=[labelPSG_all;labelPSG_seg];
        featurePSG_all=[featurePSG_all;featurePSG_seg];
        labelPSG2_all=[labelPSG2_all;labelPSG2_seg];
        featurePSG2_all=[featurePSG2_all;featurePSG2_seg];
        
        
        caseNumNCS=[caseNumNCS; ones(length(labelNCS_seg),1).*i];
        caseNumPSG=[caseNumPSG; ones(length(labelPSG_seg),1).*i];
        caseNumPSG2=[caseNumPSG2; ones(length(labelPSG2_seg),1).*i];
        
    end
    
    num=num+1;
    
end


labelNCS=labelNCS_all;
labelPSG=labelPSG_all;
labelPSG2=labelPSG2_all;

%featureNCS=permute(featureNCS_all,[1,3,2]);
featureNCS=featureNCS_all;
featurePSG=featurePSG_all;
featurePSG2=featurePSG2_all;


% [idxNan,~]=find(featurePSG~=featurePSG);
% featurePSG(idxNan,:)=[];
% labelPSG(idxNan)=[];
% caseNumPSG(idxNan)=[];    
% 
% [idxNan,~]=find(labelPSG~=labelPSG);
% featurePSG(idxNan,:)=[];
% labelPSG(idxNan)=[];
% caseNumPSG(idxNan)=[];   

%all disorder =1  all normal =0

 %LABEL 1 :snore  score=5
   %LABEL 2 :RERA score=7 + arousal score=6   
  %LABEL 3 :hyponea  score=9 + O2 desaturation score=8
 %LABEL 4 :OSA  score=10 
  %LABEL 5 :MA  score=11 
  %LABEL 6 :CSA  score=12   
  

m=[0,1,2,3,4,5,6];
m1=[0,1,1,1,1,1,1];
m2=[0,1,0,2,2,2,2];
m3=[0,0,0,1,1,1,1];
m4=[0,0,0,1,2,2,3];

labelNCS_1=relabel(labelNCS,m1);
labelPSG_1=relabel(labelPSG,m1);
labelPSG2_1=relabel(labelPSG2,m1);

labelNCS_2=relabel(labelNCS,m2);
labelPSG_2=relabel(labelPSG,m2);
labelPSG2_2=relabel(labelPSG2,m2);

labelNCS_3=relabel(labelNCS,m3);
labelPSG_3=relabel(labelPSG,m3);
labelPSG2_3=relabel(labelPSG2,m3);

labelNCS_4=relabel(labelNCS,m4);
labelPSG_4=relabel(labelPSG,m4);
labelPSG2_4=relabel(labelPSG2,m4);







Epoch_Select_ratio=sum(NCSuniqEpoch)/length(caseNumPSG);


label_count_labelNCS =uniqueCount(labelNCS);
label_count_labelPSG =uniqueCount(labelPSG);
label_count_labelPSG2 =uniqueCount(labelPSG2);

label_count_labelNCS_3 =uniqueCount(labelNCS_3);
label_count_labelPSG_3 =uniqueCount(labelPSG_3);
label_count_labelPSG2_3 =uniqueCount(labelPSG2_3);


case_count_NCS=uniqueCount(caseNumNCS);
case_count_PSG=uniqueCount(caseNumPSG);
case_count_PSG2=uniqueCount(caseNumPSG2);
ApneaRateNCS=1-(label_count_labelNCS(1,2)/sum(label_count_labelNCS(:,2)));
ApneaRatePSG=1-(label_count_labelPSG(1,2)/sum(label_count_labelPSG(:,2)));
ApneaRatePSG2=1-(label_count_labelPSG2(1,2)/sum(label_count_labelPSG2(:,2)));

ApneaRateNCS_l3=1-(label_count_labelNCS_3(1,2)/sum(label_count_labelNCS_3(:,2)));
ApneaRatePSG2_l3=1-(label_count_labelPSG2_3(1,2)/sum(label_count_labelPSG2_3(:,2)));
save([SavePath,SaveName,'.mat'],'caseNumNCS','caseNumPSG','caseNumPSG2','labelNCS','featureNCS','labelPSG','featurePSG',...
 'labelPSG2','featurePSG2','labelPSG2_2','labelPSG2_3','labelPSG2_4',...
 'labelNCS_1','labelPSG2_1','labelPSG_1',...
 'labelNCS_2','labelPSG_2','labelNCS_3','labelPSG_3','labelNCS_4','labelPSG_4','NCSuniqEpoch',...
 'label_count_labelNCS','label_count_labelPSG','label_count_labelPSG2',...
 'label_count_labelNCS_3','label_count_labelPSG_3','label_count_labelPSG2_3',...
  'case_count_NCS','case_count_PSG','case_count_PSG2',...
  'ApneaRateNCS','ApneaRatePSG','ApneaRatePSG2','trueUniq','SaveProp');


%% plot case number 
sz=13;
h1=figure();
bar(case_count_NCS(:,2));
ttl='NCS Detection Dataset';

xlabel('subject','FontSize',sz);
ylabel(['epoch number'],'FontSize',sz);
title(ttl);
set(gcf,'Position',[400,400,500,200]);
set(gca, 'FontName', 'Times New Roman');

figName1 = ['C:\Users\zz587-admin\Documents\sleep center paper\figure\matfigure\',ttl];
% print(h1,[figName1,FolderName,'.tiff'],'-dtiff','-r300');
% savefig(h1,[figName1,FolderName,'.fig']);


h1=figure();
bar(case_count_PSG2(:,2));
ttl='PSG Detection Dataset';

xlabel('subject','FontSize',sz);
ylabel(['epoch number'],'FontSize',sz);
title(ttl);
set(gcf,'Position',[400,400,500,200]);
set(gca, 'FontName', 'Times New Roman');

figName1 = ['C:\Users\zz587-admin\Documents\sleep center paper\figure\matfigure\',ttl];
% print(h1,[figName1,FolderName,'.tiff'],'-dtiff','-r300');
% savefig(h1,[figName1,FolderName,'.fig']);




function labelNew=relabel(label,m2)
labelNew=label;
m1=[0,1,2,3,4,5,6];

for i = 1: length(m1)
    labelNew(find(label==m1(i)))=m2(i);
end


end

