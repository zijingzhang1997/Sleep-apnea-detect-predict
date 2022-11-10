function [pk,VarFeature] = brEst_Pk(ncsData,fs,opts)

if ~isfield(opts,'tWinBR')
    opts.tWinBR = 8;
    fprintf('Default BR estimation window: %3.2f\n',opts.tWinBR);
end
if ~isfield(opts,'calibPk')
    opts.calibPk = 0; % Do not calibrate by default
else
    if ~isfield(opts,'calibT')
        opts.calibT = [20,40]; 
        fprintf('Default BR estimation peak calibration window [%d, %d].\n',opts.calibT(1),opts.calibT(2));
    end
    if ~isfield(opts,'calibMinPkRatio')
        opts.calibMinPkRatio = 0.4;
        fprintf('Default min peak height can be %d%% of avg peak height in calibration window.\n',100*opts.calibMinPkRatio);
    end
end

t = (0:(length(ncsData)-1))/fs;

% -------------------------------------------------------------------------
% Minima and maxima detection on NCS thorax and abdomen data
% -------------------------------------------------------------------------
pk = findMaxMin(ncsData,fs,opts);

if pk(1).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(1).ind = pk(1).ind(2:end); 
    pk(1).idx = pk(1).idx(2:end);
end

if pk(1).ind(end) == 1
    pk(1).ind = pk(1).ind(1:end-1);
    pk(1).idx = pk(1).idx(1:end-1);
end


pkMax1 = pk(1).idx(pk(1).ind == 1);
pkMin1 = pk(1).idx(pk(1).ind == 0);
par1=zeros(length(pkMax1),4);
for i=1:length(pkMax1)
    par1(i,1)=1/(t(pkMin1(i+1))-t(pkMin1(i)))*60;  %BR (BPM)
    par1(i,2)=1/2*(ncsData(pkMax1(i),1)-ncsData(pkMin1(i),1)-ncsData(pkMin1(i+1),1)); %pp average of up & down
    par1(i,3)=t(pkMax1(i))-t(pkMin1(i));  % inhale time 
    par1(i,4)=t(pkMin1(i+1))-t(pkMax1(i));  % exhale time
end

a=3; %calibrate out outlier
temp1=0;
in=zeros(length(pkMax1),1);

for i=1:length(pkMax1)    
 if abs(par1(i,1)-mean(par1(:,1))) >= a*std(par1(:,1)) ||  abs(par1(i,2)-mean(par1(:,2))) >= a*std(par1(:,2))
      
      temp1=temp1+1;
      in(i)=1;
    end
end
par1(find(in==1),:)=[];


%%
if size(ncsData,2) == 3
   if pk(2).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(2).ind = pk(2).ind(2:end); 
    pk(2).idx = pk(2).idx(2:end);
end

if pk(2).ind(end) == 1
    pk(2).ind = pk(2).ind(1:end-1);
    pk(2).idx = pk(2).idx(1:end-1);
end


pkMax2 = pk(2).idx(pk(2).ind == 1);
pkMin2 = pk(2).idx(pk(2).ind == 0);
par2=zeros(length(pkMax2),4);
for i=1:length(pkMax2)
    par2(i,1)=1/(t(pkMin2(i+1))-t(pkMin2(i)))*60;  %BR (BPM)
    par2(i,2)=1/2*(ncsData(pkMax2(i),2)-ncsData(pkMin2(i))-ncsData(pkMin2(i+1),2)); %pp average of up & down
    par2(i,3)=t(pkMax2(i))-t(pkMin2(i));  % inhale time 
    par2(i,4)=t(pkMin2(i+1))-t(pkMax2(i));  % exhale time
end

%calibrate out outlier
temp2=0;
in=zeros(length(pkMax2),1);
for i=1:length(pkMax2)
   
    if abs(par2(i,1)-mean(par2(:,1))) >= a*std(par2(:,1)) ||  abs(par2(i,2)-mean(par2(:,2))) >= a*std(par2(:,2))
      
      temp2=temp2+1;
      in(i)=1;
    end
end
par2(find(in==1),:)=[];




 if pk(3).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(3).ind = pk(3).ind(2:end); 
    pk(3).idx = pk(3).idx(2:end);
end

if pk(3).ind(end) == 1
    pk(3).ind = pk(3).ind(1:end-1);
    pk(3).idx = pk(3).idx(1:end-1);
end


pkMax3 = pk(3).idx(pk(3).ind == 1);
pkMin3 = pk(3).idx(pk(3).ind == 0);
par3=zeros(length(pkMax3),4);
for i=1:length(pkMax3)
    par3(i,1)=1/(t(pkMin3(i+1))-t(pkMin3(i)))*60;  %BR (BPM)
    par3(i,2)=1/2*(ncsData(pkMax3(i),3)-ncsData(pkMin3(i))-ncsData(pkMin3(i+1),3)); %pp average of up & down
    par3(i,3)=t(pkMax3(i))-t(pkMin3(i));  % inhale time 
    par3(i,4)=t(pkMin3(i+1))-t(pkMax3(i));  % exhale time
end

%calibrate out outlier
temp3=0;
in=zeros(length(pkMax3),1);
for i=1:length(pkMax3)
  
    if abs(par3(i,1)-mean(par3(:,1))) >= a*std(par3(:,1)) ||  abs(par3(i,2)-mean(par3(:,2))) >= a*std(par3(:,2))
      in(i)=1;
      temp3=temp3+1;
    end
end
par3(find(in==1),:)=[];

%%
  figure();
    ax(2) = subplot(3,1,1);
    plot(t,ncsData(:,1)); hold on;
    plot(t(pkMax1),ncsData(pkMax1,1),'^',...
        t(pkMin1),ncsData(pkMin1,1),'v');
    leg = {'Th','Max','Min'};
    %plotCute1('Time (s)','Resp (mV)',ax(1),[],leg,1,'Horizontal');
    ax(2) = subplot(3,1,2);
    plot(t,ncsData(:,2)); hold on;
    plot(t(pkMax2),ncsData(pkMax2,2),'^',...
        t(pkMin2),ncsData(pkMin2,2),'v');
    leg = {'Abd','Max','Min'};
    %plotCute1('Time (s)','Resp (mV)',ax(2),[],leg,1,'Horizontal');
     ax(3) = subplot(3,1,3);
    plot(t,ncsData(:,3)); hold on;
    plot(t(pkMax3),ncsData(pkMax3,3),'^',...
        t(pkMin3),ncsData(pkMin3,3),'v');
    leg = {'Abd','Max','Min'};
    

feat_mean=[mean(par1); mean(par2) ;mean(par3)];

feat_var=[var(par1); var(par2) ;var(par3)];
feat_var=feat_var./(feat_mean.^2); %normalized var ==cov  not multiple by 100 

feat_var_mean=mean(feat_var);
feat_var_var=var(feat_var);

VarFeature.mean=feat_mean;
VarFeature.var=feat_var;

BrVar=feat_var_mean(1);ppVar=feat_var_mean(2);inhaleVar=feat_var_mean(3);exhaleVar=feat_var_mean(4);
table1=table(BrVar,ppVar,inhaleVar,exhaleVar);
VarFeature.var_mean=table1;
BrVar=feat_var_var(1);ppVar=feat_var_var(2);inhaleVar=feat_var_var(3);exhaleVar=feat_var_var(4);
table2=table(BrVar,ppVar,inhaleVar,exhaleVar);
VarFeature.var_var=table2;

fprintf('calibrate out %d %d %d breath cycle\n',temp1,temp2,temp3);
else
%% Only one NCS column    
    % Find breath rate in last tWinBR sec (or slightly less)
    br = zeros(length(t),1);
    var1=zeros(length(t),1);  % BR var
    var2=zeros(length(t),1);  %pp var
    var3=zeros(length(t),2);  %inhale exhale time var
  
    t = t(:);
    tBRmax1 = t(pkMax2);

    for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinBR))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            br(iter,1) = 0;
        else  
            %average calculatio of BR
            for i=1:length(idxBR)-1
               br(iter,1) = br(iter,1)+(idxBR(i+1)- idxBR(i))/(tBRmax1(idxBR(i+1))-tBRmax1(idxBR(i)))/(length(idxBR)-1);   
            end       
            % br(iter,1) = (idxBR(2)-idxBR(1))/(tBRmax1(idxBR(2))-tBRmax1(idxBR(1)));
            %only use one cycle    
        end
    end
    %BR variation
    for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinVar))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            var1(iter,1) = 0;
        else  
            %average calculatio of BR
            temp=zeros(length(idxBR)-1,1);
            
            for i=1:length(idxBR)-1
               temp(i) = (idxBR(i+1)- idxBR(i))/(tBRmax1(idxBR(i+1))-tBRmax1(idxBR(i))); 
               temp(i)=temp(i)*60;
            end
            var1(iter,1)=var(temp)/((mean(temp)).^2);
            mean1(iter,1)=mean(temp);
            % br(iter,1) = (idxBR(2)-idxBR(1))/(tBRmax1(idxBR(2))-tBRmax1(idxBR(1)));
            %only use one cycle    
        end
    end 
    %peak-peak variation
       for iter = 1:length(t)

        idxBR = find((tBRmax1 > (t(iter)-opts.tWinpp))&(tBRmax1 <= t(iter)));

        if (size(idxBR,1) <= 1)
            var2(iter,1) = 0;
        else  
            %average calculatio of BR
            temp=zeros(2*(length(idxBR)-1),1);
            temp_in=zeros(length(idxBR)-1,1);
            temp_ex=zeros(length(idxBR)-1,1);
            for i=1:length(idxBR)-1
                a1=pkMax2(idxBR(i));
                a2=pkMax2(idxBR(i+1));
                b=pkMin2(idxBR(i+1));
                temp(2*i-1) = ncsData(a2)-ncsData(b);
                temp_in(i)=t(a2)-t(b);
                temp_ex(i)=t(b)-t(a1);
                temp(2*i) = ncsData(a1)-ncsData(b); 
               
            end
            mean2(iter,1)=mean(temp);mean3(iter,1)=mean(temp_in);mean4(iter,1)=mean(temp_ex);
            var2(iter,1)=var(temp)/((mean(temp))^2);%CoV = std/mean
            var3(iter,1)=var(temp_in)/((mean3(iter,1))^2);  %normalize??
            var3(iter,2)=var(temp_ex)/((mean4(iter,1))^2);
           
        end
       end 
 var2=var2.*100; 
 var1=var1.*100; 
 var3=var3.*100; 
 Var=[var1,var2,var3];
    br = 60.*br;
    
meanBR=mean(mean1);meanPP=mean(mean2);meanIn=mean(mean3);meanEx=mean(mean4);
VarBR_mean=mean(var1(var1~=0));
VarBR_std=std(var1(var1~=0));
VarPP_mean=mean(var2(var2~=0));
VarPP_std=std(var2(var2~=0));
VarIn=var3(:,1);
VarIn_mean=mean(VarIn(VarIn~=0));
VarIn_std=std(VarIn(VarIn~=0));
VarEx=var3(:,2);
VarEx_mean=mean(VarEx(VarEx~=0));
VarEx_std=std(VarEx(VarEx~=0));

VarFeature.VarBR_mean=VarBR_mean;
VarFeature.VarBR_std=VarBR_std;
VarFeature.Varpp_mean=VarPP_mean;
VarFeature.Varpp_std=VarPP_std;
VarFeature.VarIn_mean=VarIn_mean;
VarFeature.VarIn_std=VarIn_std;
VarFeature.VarIn_mean=VarIn_mean;
VarFeature.VarIn_std=VarIn_std;
VarFeature.VarEx_mean=VarEx_mean;
VarFeature.VarEx_std=VarEx_std;
VarFeature.meanBR=meanBR;VarFeature.meanPP=meanPP;VarFeature.meanIn=meanIn;VarFeature.meanEx=meanEx;

    figure();
    ax(1) = subplot(1,1,1);
    plot(t,ncsData); hold on;
    plot(t(pkMax2),ncsData(pkMax2,1),'^',...
        t(pkMin2),ncsData(pkMin2,1),'v');
    leg = {'Resp','Max','Min'};
  plotCute1('Time (s)','Resp (mV)',ax(1),[],leg,1,'Horizontal');
%     ax(2) = subplot(nFig,1,2);
%     plot(t,br(:,1)); 
%     leg = {'BR'};
%     plotCute1('Time (s)','BR (BPM)',ax(2),[],leg,1,'Horizontal');
%     linkaxes(ax,'x');
%     
end

