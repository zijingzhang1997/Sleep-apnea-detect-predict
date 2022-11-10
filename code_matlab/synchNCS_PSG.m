function [NCSdata,PSGdata,score,score_all,StTime]=synchNCS_PSG(NCSdata,StTimeNCS,PSGdata,StTimePSG,scorePSG,scorePSG_all,fs)


score=scorePSG;
score_all=scorePSG_all;
EndTimeNCS=StTimeNCS +seconds(size(NCSdata,1)/fs);
EndTimePSG=StTimePSG +seconds(size(PSGdata,1)/fs);

if EndTimeNCS<StTimePSG | EndTimePSG<StTimeNCS
   fprintf('no mutual time \n')
   return
end
StTime=max(StTimeNCS,StTimePSG);
EndTime=min(EndTimeNCS,EndTimePSG);

if StTime>StTimeNCS 
    NCSdata(1:seconds(StTime-StTimeNCS)*fs,:)=[];
    
end

if  StTime>StTimePSG
    PSGdata(1:seconds(StTime-StTimePSG)*fs,:)=[];
    score(1:seconds(StTime-StTimePSG)*fs,:)=[];
    score_all(1:seconds(StTime-StTimePSG)*fs,:)=[];
end


if EndTime<EndTimeNCS 
    del=-seconds(EndTime-EndTimeNCS);
    NCSdata(end-del*fs+1:end,:)=[];
    
end
if  EndTime<EndTimePSG
    del=-seconds(EndTime-EndTimePSG);
    PSGdata(end-del*fs+1:end,:)=[];
    score(end-del*fs+1:end,:)=[];
    score_all(end-del*fs+1:end,:)=[];
end   
end
