function [ch1234Ds,StartTime]=readNCSmat(filePathName,fsDS)
load(filePathName);
% The order is [amp_tx1rx1 amp_tx2rx2 ph_tx1rx1 ph_tx2rx2]

ch1234 = [ ConvertedData.Data.MeasuredData(4).Data,...
        ConvertedData.Data.MeasuredData(5).Data,...
        ConvertedData.Data.MeasuredData(6).Data,...
        ConvertedData.Data.MeasuredData(7).Data,... 
        ];
StartTime=ConvertedData.Data.MeasuredData(4).Property(1).Value; % start time in datetime format.
StartTime=datetime(StartTime);

%StarTime = datevec(StartTime);

ch1234Ds=resample(ch1234,fsDS,500,'Dimension',1);


end