dataPath1=['C:\Users\zz587-admin\Documents\sleep center\data_reference\NCS EDF records\'];


for i=24:25
fileName=['NCS0',num2str(i,'%02d')];
filePathName= [dataPath1,fileName,'.edf'];

info = edfinfo(filePathName);

startTime=append(info.StartDate,' ',info.StartTime);
StTime(i) = datetime(startTime,'InputFormat','dd.MM.yy HH.mm.ss');

end


StTimeVec=datevec(StTime);