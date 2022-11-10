function [waveDS,waveSpO2DS,StTimeRec,score,score_all,annotations]=readPSG(filePathName,fsDS)
%waveDS  downsampled to fsDS waveform *3 channels ("Airflow" "Chest" "Abdomen")  len=samplePoints=sec*fsDS
%StTimeRec= startTime (waveform start)+ AddSecond = score start time.  datetime format
%score =[time  label]  len=timelength(s) label 0-10

fs=250;
[dataall,annotations] = edfread(filePathName);
info = edfinfo(filePathName);

startTime=append(info.StartDate,' ',info.StartTime);
StTime = datetime(startTime,'InputFormat','dd.MM.yy HH.mm.ss');
%StTimeVec = datevec(StTime);

[data,annotations] = edfread(filePathName);
if contains(join(info.SignalLabels),"Airflow")
waveTT=edfread(filePathName,'SelectedSignals',["Airflow" "Chest" "Abdomen"]);
else 
 fprintf('no Airflow Channel ')   
waveTT=edfread(filePathName,'SelectedSignals',[ "Chest" "Abdomen"]);


end
waveSpO2=edfread(filePathName,'SelectedSignals',[ "SpO2"]);


 
%get waveforms (for calibreation of NCS)and downsample
TTable = timetable2table(waveTT);
T=table2array(TTable(:,2:end));
T=cell2mat(T);
waveDS=resample(T,fsDS,fs,'Dimension',1);

% if no air flow substitue with chest 
if ~contains(join(info.SignalLabels),"Airflow")
    fprintf('substitute airflow with chest ') 
    waveDS=[waveDS(:,1) waveDS(:,2) waveDS(:,1)];
end

%get other waveforms spo2
SpO2Table=timetable2table(waveSpO2);
SpO2=cell2mat(table2array(SpO2Table(:,2:end)));
SpO2(SpO2==0)=mean(SpO2);
waveSpO2DS=resample(SpO2,fsDS,25,'Dimension',1);      %fs sp02=25



StAddSec=waveTT.Properties.StartTime;
StTimeRec=StTime+StAddSec;  % start time in datetime format.  starttime+startSecond in recording 


   
TTime=[0:size(TTable,1)-1];  % record time array (in sec) -startRecordTime
TTime=TTime';
t=linspace(1/fsDS,size(TTable,1),fsDS*size(TTable,1));
[score,score_all]=readScore(annotations,t,fsDS);  %change: score time unit(sec/fs)

end

