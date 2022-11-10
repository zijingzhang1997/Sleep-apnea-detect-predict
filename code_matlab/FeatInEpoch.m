% feature=[mean(br), std(br), mean(pp), std(pp), mean(in), std(in), mean(ex), std(ex),... 
%          skew_mean, kurt_mean, entro,...
%            per_power, cycle, covBR, covPP, void_t, ...
%            max_in, max_ex, max_br, min_br, min_pp, ...
%            Cor_br, SD_br, Cor_pp, SD_pp, Cor_in, SD_in, Cor_ex, SD_ex];



function featAll=FeatInEpoch(EpochData,fs,opt)
featAll=zeros(opt.featNum,size(EpochData,2)); 

for i=1:size(EpochData,2)
    EpochData_row=EpochData(:,i);
    
    featAll(:,i)=featExtract(EpochData_row,fs,opt);
    
end  

end

    

function feature=featExtract(Data,fs,opt)    
    

%feature=zeros(opt.featNum,1);   
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
f(1:5)=[];pxx(1:5)=[];
fmax=f(find(pxx==max(pxx)));
covBR=std(br)/mean(br);
covPP=std(pp)/mean(pp);

max_br=max(br);
min_br=min(br);
max_in=max(in);
max_ex=max(ex);
min_pp=min(pp);
if isempty(br)
    max_br=0;
    min_br=0;
end
if isempty(pp)
    min_pp=0;
    
end
if isempty(in) || isempty(ex) 
    max_in=0;
    max_ex=0;
end


% auto correlation and successive difference 
Cor_br = xcorr(br,1,'coeff');Cor_br =Cor_br (1);
SD_br= mean(abs(diff(br))./br(1:end-1));
Cor_pp = xcorr(pp,1,'coeff');Cor_pp =Cor_pp (1);
SD_pp= mean(abs(diff(pp))./pp(1:end-1));
Cor_in = xcorr(in,1,'coeff');Cor_in =Cor_in (1);
SD_in= mean(abs(diff(in))./in(1:end-1));
Cor_ex = xcorr(ex,1,'coeff');Cor_ex =Cor_ex (1);
SD_ex= mean(abs(diff(ex))./ex(1:end-1));

% feature=[mean(br) std(br) mean(pp) std(pp) mean(in) std(in) mean(ex) std(ex)... 
%          skew_mean kurt_mean entro...
%            per_power cycle covBR covPP void_t ];
     
feature=[mean(br), std(br), mean(pp), std(pp), mean(in), std(in), mean(ex), std(ex),... 
         skew_mean, kurt_mean, entro,...
           per_power, cycle, covBR, covPP, void_t, ...
           max_in, max_ex, max_br, min_br, min_pp, ...
           Cor_br, SD_br, Cor_pp, SD_pp, Cor_in, SD_in, Cor_ex, SD_ex];
     
feature(isnan(feature))=0;




    
end



