function [featAll,pkidx]=PlotMaxMin_FeatInEpoch(EpochData,fs,opt)

for i=1:size(EpochData,2)
    EpochData_row=EpochData(:,i);
    
    [featAll(:,i),pkidx(i)]=PlotMaxMin_featExtract(EpochData_row,fs,opt);
    
        
end  

end

    
function [feature,pkidx]=PlotMaxMin_featExtract(Data,fs,opt)    
    

    
t = (0:(length(Data)-1))/fs;
Data = rescale(Data,-1,1); 
% Minima and maxima detection 
pk = findMaxMin(Data,fs,opt);

if pk(1).ind(1) == 1
    % If first peak is inhalation peak, skip it
    pk(1).ind = pk(1).ind(2:end); 
    pk(1).idx = pk(1).idx(2:end);
end

if pk(1).ind(end) == 0
    pk(1).ind = pk(1).ind(1:end-1);
    pk(1).idx = pk(1).idx(1:end-1);
end
pkMax = pk(1).idx(pk(1).ind == 1);
pkMin = pk(1).idx(pk(1).ind == 0);

%% calibrate peaks 
if opt.calib==1
  
   
    del = zeros(length(pkMax),1);

   
    for i = 1:length(pkMax)
            del(i) = abs(Data(pkMax(i))-Data(pkMin(i))); % Making it positive always
    end

    tBRmax1 = t(pkMax);
    idxCalib = ((tBRmax1 >= opt.calibT(1))&(tBRmax1 <= opt.calibT(2)));
    delNcsCalib = mean(del(idxCalib)); 


    idxPk = (del >= opt.calibMinPkRatio*delNcsCalib);
    nold=length(pkMax);
    pkMax = pkMax(idxPk); % Only keep indices that satisfy, otherwise ignore
    pkMin = pkMin(idxPk);
    nn=length(pkMax);
    CaliDel=nold-nn;
    
end

    %% Find features in each cycle (min-max-min)
cycle=length(pkMax)-1;
br = zeros(cycle,1);pp = zeros(cycle,1);in = zeros(cycle,1);ex = zeros(cycle,1);
skew=zeros(cycle,1);kurt=zeros(cycle,1);
for i=1:cycle
        br(i)=60/(t(pkMin(i+1))-t(pkMin(i)));  %BR unit BPM 
        pp(i)=(Data(pkMax(i))-Data(pkMin(i))-Data(pkMin(i+1)))/2;
        in(i)=t(pkMax(i))-t(pkMin(i));
        ex(i)=t(pkMin(i+1))-t(pkMax(i));
        
        DataSeg=Data(pkMin(i):pkMin(i+1));
        kurt(i)=kurtosis(DataSeg);  %measure of the "tailedness"
        skew(i)=skewness(DataSeg);  %measure of the asymmetry
        en(i)=entropy(DataSeg);
    end
 
entro=entropy(Data);  %entropy of waveform in entire window > entropy in each waveform cycle ?
skew_mean=mean(skew);kurt_mean=mean(kurt);
void_t=t(end)-sum(in)-sum(ex);  % void  time: duration  no peak is detected 
  
[pxx,f] = periodogram(Data,hamming(length(Data)),length(Data),fs,'power');

pBand = bandpower(pxx,f,[opt.psdBand(1) opt.psdBand(2)],'psd');
pTot = bandpower(pxx,f,'psd');
per_power = 100*(pBand/pTot);

fmax=f(find(pxx==max(pxx)));

feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
         skew_mean kurt_mean entro...
           per_power cycle fmax CaliDel void_t];
pkidx.max=pkMax;    
pkidx.min=pkMin; 





    
end






