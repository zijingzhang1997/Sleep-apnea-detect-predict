function   [result,result_all]=readScore(annotations,TTime,fs);
% Sleep stage   N1=1 N2=2 R=4 W=3
% Snore=5 
% Arousal=6 RERA=7
% O2 Desaturation=8  Hypopnea=9
%  Obstructive Apnea=10 Mixed Apnea=11 Central Apnea=12
%   %   tcpCO2 LM Impedance ??
%result= [time  label]  timelength(s)*2    time(unit:sec) start from 0 second (startTimeSecond + dateTime ) to end.
% Label 0-11. 1-4 sleep stage N1 N2 R W ; 5-10 apnea ; 0 no score    

label=zeros(length(TTime),1);
label_all=zeros(length(TTime),12);   %back up label in case of repeat score on same time sample

stage=annotations.Annotations;
duration=annotations.Duration;
duration=seconds(duration);% duration in seconds 
RowTimes=annotations.Properties.RowTimes;
RowTimes=RowTimes(:)-annotations.Properties.StartTime;%time stamp - startTime  in duration format
RowTimes=seconds(RowTimes); % startTime  in seconds
stg=["Sleep stage N1","Sleep stage N2","Sleep stage R","Sleep stage W",...
    "Snore","Arousal","RERA","O2 Desaturation","Hypopnea","Obstructive Apnea","Mixed Apnea","Central Apnea"];
lb=[1:12];
duration(duration~=duration)=0;

for i=1:size(annotations,1)
    for j=1:length(stg)
    if contains(stage(i),stg(j))
        
        ind=[round(RowTimes(i)*fs)+1:round((RowTimes(i)+duration(i))*fs)];
        if ind~=NaN
        label(ind)=lb(j);
        label_all(ind,j)=1;
        end
        
    end
    end
end

for i=1:size(label_all,1)
    idx=find(label_all(i,:));
    if isempty(idx) 
        label(i)=0;
    else
        label(i)=max(idx);
    end
end







l=min([length(TTime);length(label)]) ;   % revised due to error 
result=[TTime(1:l)' label(1:l,:)];
result_all=[TTime(1:l)' label_all(1:l,:)];

[C,ia,ic] = unique(label);
a_counts = accumarray(ic,1);
label_counts = [C, a_counts];

end