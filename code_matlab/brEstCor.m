%revise by zijing from 'brEst'
% find all peaks , calculate correaltion of all peaks with one unit of gap 
%plot the scatter plot of x(current BR), y(next BR)
function [cor,h,ncsDataNew,temp,SD] = brEstCor(ncsData,fs,CaseNum)

opts.calibPk = 1;
if ~isfield(opts,'calibPk')
    opts.calibPk = 0; % Do not calibrate by default
else
    if ~isfield(opts,'calibT')
        opts.calibT = [20,40]; 
      %  fprintf('Default BR estimation peak calibration window [%d, %d].\n',opts.calibT(1),opts.calibT(2));
    end
    if ~isfield(opts,'calibMinPkRatio')
        opts.calibMinPkRatio = 0.4;
      %  fprintf('Default min peak height can be %d%% of avg peak height in calibration window.\n',100*opts.calibMinPkRatio);
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


pkMax = pk(1).idx(pk(1).ind == 1);
pkMin = pk(1).idx(pk(1).ind == 0);

if opts.calibPk == 1
   % fprintf('\nBR: Performing calibration of pk-pk height.\n')
    % This is the change in Pk-Pk thorax and abdomen signal
    del2 = zeros(length(pkMax),1);

    % So for a cycle, considering there exists 2 minima and 2 maxima point:
    % Calculation is peformed using the difference between maxima and first
    % minima. The update is performed at the end of cycle to be consistent with
    % TV from airflow calculation. And this value is held until next update or
    % the end of the waveform.
    for i = 1:length(pkMax)
            del2(i) = abs(ncsData(pkMax(i),1)-ncsData(pkMin(i),1)); % Making it positive always
    end

    tBRmax1 = t(pkMax);
    idx2Calib = ((tBRmax1 >= opts.calibT(1))&(tBRmax1 <= opts.calibT(2)));
    delNcs2Calib = mean(del2(idx2Calib)); 
%     fprintf('delNCS2Calib = %f\n',delNcs2Calib);

    idxPk = (del2 >= opts.calibMinPkRatio*delNcs2Calib);

    pkMax = pkMax(idxPk); % Only keep indices that satisfy, otherwise ignore
    pkMin = pkMin(idxPk);
    
    pk(1).idxValidPk = idxPk;
end

% fig1 = figure;  by zijing
% nFig = size(ncsData,2)+1;

%%



%% Only one NCS column    
    % Find breath rate in last tWinBR sec (or slightly less)
    br = zeros(length(pkMin)-1,1);
  
  
    t = t(:);
    tmax = t(pkMax);
    tmin = t(pkMin);
    br=(1./diff(tmin))*60;  
    for i=1:length(pkMax)-1
        pp(2*i-1)=ncsData(pkMax(i))-ncsData(pkMin(i));
        in(i)=t(pkMax(i))-t(pkMin(i));
        pp(2*i)=ncsData(pkMax(i))-ncsData(pkMin(i+1));
        ex(i)=t(pkMin(i+1))-t(pkMax(i));
       
    end
    for i=length(pkMax)   
        pp(2*i-1)=ncsData(pkMax(i))-ncsData(pkMin(i));
        in(i)=t(pkMax(i))-t(pkMin(i));
    end
   %% delete outlier
  

switch CaseNum(2)
    case 3
        a=3.5;
    case 5
        a=3.5;
    case 1
        a=3;    
    otherwise
        a=3;
end

        
        
        temp=0;
        ind=zeros(length(pkMax),1);
      for i=1:length(pkMax)-1
          if abs(br(i)-mean(br)) >= a*std(br) || abs(pp(2*i-1)-mean(pp)) >= a*std(pp) || abs(pp(2*i)-mean(pp)) >= a*std(pp)...
                  || abs(in(i)-mean(in)) >= a*std(in)|| abs(ex(i)-mean(ex)) >= a*std(ex)
              ind(i)=1;
              temp=temp+1;
          end
      end
      tInd=find(ind==1);
      tInd_del=[];
      for j=1:temp
         t_temp= [pkMin(tInd(j)):pkMin(tInd(j)+1)];
         tInd_del=[tInd_del  t_temp];
      end
      ncsData(tInd_del)=[];
      br(tInd)=[];in(tInd)=[];ex(tInd)=[];
      ppInd=[tInd.*2;tInd.*2-1];pp(ppInd)=[];
      ncsDataNew=ncsData;
      fprintf('calibrate out %d  breath cycle\n',temp);
   
SD=[std(br);std(pp);std(in);std(ex);mean(br);mean(pp);mean(in);mean(ex)];

SD=round(SD,4);
      
  %% calculate cor  
    brCor1 = xcorr(br,1,'coeff');  brCor2= autocorr(br,'NumLags',1);brCor3= mean(abs(diff(br))./br(1:end-1));brCor=[brCor1(1) brCor2(2) brCor3];
    ppCor1 = xcorr(pp,1,'coeff');  ppCor2= autocorr(pp,'NumLags',1);ppCor3= mean(abs(diff(pp))./pp(1:end-1));ppCor=[ppCor1(1) ppCor2(2) ppCor3];
    inCor1 = xcorr(in,1,'coeff');  inCor2= autocorr(in,'NumLags',1);inCor3= mean(abs(diff(in))./in(1:end-1));inCor=[inCor1(1) inCor2(2) inCor3];
    exCor1 = xcorr(ex,1,'coeff');  exCor2= autocorr(ex,'NumLags',1);exCor3= mean(abs(diff(ex))./ex(1:end-1));exCor=[exCor1(1) exCor2(2) exCor3];
    cor=[brCor; ppCor; inCor; exCor];
    cor=round(cor,4);
 sz=10;
h=figure;
subplot(2,2,1)
scatter(br(1:end-1),br(2:end),'*','r');
xlim([0 50]);ylim([0 50]);
title(['r_1=',num2str(cor(1,1)),' r_2=',num2str(cor(1,3)),' r_3=',num2str(cor(1,3))],'FontSize',sz)
xlabel('current BR (BPM)','FontSize',sz);
ylabel('next BR (BPM)','FontSize',sz);
a=0:5:60;
hold on
plot(a,a,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);

subplot(2,2,2)
scatter(pp(1:end-1),pp(2:end),'*','b');
xlim([0 2]);ylim([0 2]);
title(['r_1=',num2str(cor(2,1)),' r_2=',num2str(cor(2,3))],'FontSize',sz)
xlabel('current peak-peak value (a.u.)','FontSize',sz);
ylabel('next peak-peak value (a.u.)','FontSize',sz);
a=0:0.5:5;
hold on
plot(a,a,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);

subplot(2,2,3)
scatter(in(1:end-1),in(2:end),'*','g');
xlim([0 4]);ylim([0 4]);
title(['r_1=',num2str(cor(3,1)),'  r_2=',num2str(cor(3,3))],'FontSize',sz)
xlabel('current inhalation time (s)','FontSize',sz);
ylabel('next inhalation time (s)','FontSize',sz);
a=0:0.5:5;
hold on
plot(a,a,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);

    subplot(2,2,4)
scatter(ex(1:end-1),ex(2:end),'*','c');
xlim([0 4]);ylim([0 4]);
title(['r_1=',num2str(cor(4,1)),'  r_2=',num2str(cor(4,3))],'FontSize',sz)
xlabel('current exhalation time (s)','FontSize',sz);
ylabel('next exhalation time (s)','FontSize',sz);
a=0:0.5:5;
hold on
plot(a,a,'color',[0.5,0.5,0.5],'LineStyle',':','LineWidth',2);

RoutineName=['Case',num2str(CaseNum(1)),'Routine',num2str(CaseNum(2))];
sgtitle(RoutineName)
set(gcf,'Position',[100,200,500,480]);
    
%     figure();
%     ax(1) = subplot(1,1,1);
%     plot(t,ncsData); hold on;
%     plot(t(pkMax),ncsData(pkMax,1),'^',...
%         t(pkMin),ncsData(pkMin,1),'v');
%     leg = {'Resp','Max','Min'};
%   plotCute1('Time (s)','Resp (mV)',ax(1),[],leg,1,'Horizontal');
  

end

