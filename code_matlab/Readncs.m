dataPath=['C:\Users\zz587-admin\Documents\sleep center\matFile\sc1\'];
fileName='sc1_ 1012_0934'
fs=250e3;
fsDS=500;
toff=[10:230]';

filePathName = [dataPath,fileName,'.tdms'];
filePathName_m = [dataPath,fileName,'.mat'];
if ~exist(filePathName_m,'file') & ~exist(filePathName,'file')
    return
end
if ~exist(filePathName_m,'file')
   convertTDMS_sim(true,filePathName);
end

load(filePathName_m);

ampTh=ConvertedData.Data.MeasuredData(4).Data;
ampAb=ConvertedData.Data.MeasuredData(5).Data;
phTh=ConvertedData.Data.MeasuredData(6).Data;
phAb=ConvertedData.Data.MeasuredData(7).Data;

ampdsTh=resample(ampTh,fsDS,fs);
ampdsAb=resample(ampAb,fsDS,fs);
t = ((0:(length(ampdsTh)-1))/fsDS)';
ampfiltTh=ampdsTh((toff(1)*fsDS):toff(size(toff))*fsDS);
ampfiltAb=ampdsAb((toff(1)*fsDS):toff(size(toff))*fsDS);
tOff=((0:(length(ampfiltTh)-1))/fsDS)';
opts1.filtType = 'LpHp'; opts1.orderHP = 5;
opts1.f3db = 0.05; opts1.fpLP = 1; opts1.fstLP = 1.5;
ampfiltTh = filterLpHp(ampfiltTh,fsDS,opts1); % th amp
ampfiltAb = filterLpHp(ampfiltAb,fsDS,opts1); 



tOff=((0:(length(ampfiltTh)-1))/fsDS)';

% opts1.filtType = 'LpHp'; opts1.orderHP = 5;
% opts1.f3db = 0.05; opts1.fpLP = 5; opts1.fstLP = 7;





[ampfiltAb,PS] = mapminmax(ampfiltAb');
ampfiltAb=ampfiltAb';
[ampfiltTh,PS] = mapminmax(ampfiltTh');
ampfiltTh=ampfiltTh';


figure
sz=10;
subplot(2,1,1);
plot(tOff,ampfiltTh,'LineWidth',0.5,'color','red');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])

title('Th sitting 0.9GHz','FontSize',sz)

subplot(2,1,2);
plot(tOff,ampfiltAb,'LineWidth',0.5,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('Amp (a.u.)','FontSize',sz)
xlim([0 220])


title('Ab sitting 0.9GHz','FontSize',sz)










%%  breath rate estimation  peak-to-peak
opts3.tWinBR = 15; % Window on which br is estimated
opts3.tWin = 4; % Window for peak detection moving average
opts3.minInterceptDist = 0.15; 
opts3.calibPk = 1; % Calibrate out peaks of small height (max-min)
opts3.calibMinPkRatio = 0.5;
ampfiltTh=cat(1,ampfiltTh(1:toff2(1)*fsDS),ampfiltTh(toff2(2)*fsDS:length(ampfiltTh)));
tOff=((0:(length(ampfiltTh)-1))/fsDS)';
ampfiltAb=cat(1,ampfiltAb(1:toff2(1)*fsDS),ampfiltAb(toff2(2)*fsDS:length(ampfiltAb)));
[brTh,ThRespPk] = brEst(ampfiltTh,fsDS,opts3);
pkMaxTh = ThRespPk.idx(ThRespPk.ind == 1);
pkMinTh = ThRespPk.idx(ThRespPk.ind == 0);
pkMaxTh = pkMaxTh(ThRespPk.idxValidPk);
pkMinTh = pkMinTh(ThRespPk.idxValidPk);

[brAb,AbRespPk] = brEst(ampfiltAb,fsDSbio,opts3);
pkMaxAb = AbRespPk.idx(AbRespPk.ind == 1);
pkMinAb = AbRespPk.idx(AbRespPk.ind == 0);
pkMaxAb = pkMaxAb(AbRespPk.idxValidPk);
pkMinAb = pkMinAb(AbRespPk.idxValidPk);


% rmseBR = sqrt(mean((brBio-brNcs).^2));
% fprintf('Done. RMSE BR: %3.2f\n',rmseBR);

sz=10;
zero_long=zeros(toff2(2)*fsDS-toff2(1)*fsDS,1);
brTh_long=cat(1,brTh(1:toff2(1)*fsDS),zero_long,brTh(toff2(1)*fsDS +1:length(brTh)));
brAb_long=cat(1,brAb(1:toff2(1)*fsDS),zero_long,brAb(toff2(1)*fsDS +1:length(brAb)));
tOfflong=((0:(length(brAb_long)-1))/fsDS)';

figure()
plot(tOfflong,brTh_long,'LineWidth',1,'color','red');

hold on
plot(tOfflong,brAb_long,'LineWidth',1,'color','green');
xlabel('time (s)','FontSize',sz)
ylabel('BR (BPM)','FontSize',sz)
xlim([0 220])
legend('Th','Ab','FontSize',sz)
%%
% heart beat rate filter
% yampDsOff=fft(ampfiltHRTh);
% f = (0:length(yampDsOff)-1)*fsDS/length(yampDsOff)';
% yam=abs(yampDsOff);
% f=f(1,1:1000);
% yam=yam(1:1000,1)';
% yammax=max(yam);
% yammin=min(yam);
% yamnorm=zeros(1,1000);
% yamnorm=(yam-yammin)/(yammax-yammin);
% figure
% plot(f,yamnorm,'LineWidth',0.5,'color','blue');
% xlabel('frequency/Hz')
% 
% title('amp Th')


% % %%  Using cross correlation to estimate time shift, in case
% tCorr = [5, 50]; 
% fprintf('\nTime sync, Calibration: correlating window [%d,%d]s\n',tCorr(1),tCorr(2));
% 
% nStart = tCorr(1)*fsDS; % Assuming same fs for both ncs and biopac
% nEnd = tCorr(2)*fsDS; 
% maxLag = 1000;
% [r, lags] = xcorr(ampfiltThnorm(nStart:nEnd),ampfiltAbnorm(nStart:nEnd),maxLag);
% figure; plot(lags,r)
% [~,rMaxIdx] = max(abs(r));
% lagsMax = lags(rMaxIdx);
% tDevCalib = lagsMax/fsDS;
% fprintf('Suggested NCS Th-Bio Th calibration time offset is %f\n',tDevCalib);   
%  
% 
% 
% %% ONLY call this block once!!
% % Now shift the waveforms to compensate the time difference
% % Only shift if time difference is more than a threshold and less than a
% % threshold: Ideally minimum should be sampling frequency
% 
% fprintf('Performing synchronization based on time-shift estimate...  \n ');
% thMinCorr = 0.4;
% tOffMinMax = [0.005, 1.0];
% % First the calibration waveforms (Bio, NCS Th):
% if abs(tDevCalib)>=tOffMinMax(1) && abs(tDevCalib) <= tOffMinMax(2) 
%     nSampDev = abs(tDevCalib) * fsDS; % same BIO and NCS frequencies
%     if tDevCalib > 0
%         ampfiltThnorm = ampfiltThnorm(nSampDev+1:end,:);
%         ampfiltAbnorm = ampfiltAbnorm(1:end-nSampDev,:);
%     else
%         ampfiltThnorm = ampfiltThnorm(1:end-nSampDev,:);
%         ampfiltAbnorm = ampfiltAbnorm(nSampDev+1:end,:);
%     end
%     toffcorr = tOff(1:end-nSampDev);
%     fprintf('\nSync: NCS calib Th with Biopac.\n');
% else
%     tDevCalib = 0; % To not be recorded, since we are not correcting
% end
% 
% figure()
% plot(toffcorr,ampfiltThnorm,'LineWidth',0.5,'color','red');
% 
% hold on
% plot(toffcorr,ampfiltAbnorm,'LineWidth',0.5,'color','green');
% xlabel('t/s')
% ylabel('normalized Amplitude')
% legend('amp NCS TH','amp NCS AB')