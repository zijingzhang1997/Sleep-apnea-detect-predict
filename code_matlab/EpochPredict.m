function labelpred=EpochPredict(label,opt)
labelpred=zeros(length(label),1);
%% Not prediction, only detection. label the epoch as it is now 
if opt==0
    
    labelpred=label;
    return 
end


%% Apnea epoch included. As long as epoch after it is apnea. we label it as
% apnea 
if opt==1
    for i=1:length(label)-2
        
        if  label(i+1)>0
            labelpred(i)=label(i+1);
        elseif  label(i+2)>0
            labelpred(i)=label(i+2);
        end
        
        
        % if label itself is apnea but the epochs after is normal it will be
        % treated as normal . So use NAN
        if label(i)>0 & label(i+1)==0 & label(i+2)==0
            labelpred(i)=nan;
        end
    end
end


%% Apnea epoch excluded. Only use normal epoch to predict 
if opt==2
    for i=1:length(label)-2
        
        if  label(i+1)>0
            labelpred(i)=label(i+1);
        elseif  label(i+2)>0
            labelpred(i)=label(i+2);
        end
        
        
        % if label itself is apnea, we directly delete it 
        % . So use NAN
        if label(i)>0 
            labelpred(i)=nan;
        end
    end
end

%% Apnea epoch excluded. Only use normal epoch to predict 
% in advance of 4   epoch   consider overlap of epochs 
if opt==3
    for i=1:length(label)-4
        
        if  label(i+1)>0
            labelpred(i)=label(i+1);
        elseif  label(i+2)>0
            labelpred(i)=label(i+2);
        elseif  label(i+3)>0
            labelpred(i)=label(i+3);
        elseif  label(i+4)>0
            labelpred(i)=label(i+4);
        end
        
        
        
        % if label itself is apnea, we directly delete it 
        % . So use NAN
        if label(i)>0 
            labelpred(i)=nan;
        end
    end
end
%% Apnea epoch excluded. Only use normal epoch to predict 
% in advance of 6   epoch   consider overlap of epochs 
if opt==4
    for i=1:length(label)-6
        
        if  label(i+1)>0
            labelpred(i)=label(i+1);
        elseif  label(i+2)>0
            labelpred(i)=label(i+2);
        elseif  label(i+3)>0
            labelpred(i)=label(i+3);
        elseif  label(i+4)>0
            labelpred(i)=label(i+4);
         elseif  label(i+5)>0
            labelpred(i)=label(i+5);
         elseif  label(i+6)>0
            labelpred(i)=label(i+6);
        end
        
        
        
        % if label itself is apnea, we directly delete it 
        % . So use NAN
        if label(i)>0 
            labelpred(i)=nan;
        end
    end
end
%% Apnea epoch excluded. Only use normal epoch to predict 
% in advance of 8   epoch   consider overlap of epochs 
if opt==5
    for i=1:length(label)-8
        
        if  label(i+1)>0
            labelpred(i)=label(i+1);
        elseif  label(i+2)>0
            labelpred(i)=label(i+2);
        elseif  label(i+3)>0
            labelpred(i)=label(i+3);
        elseif  label(i+4)>0
            labelpred(i)=label(i+4);
         elseif  label(i+5)>0
            labelpred(i)=label(i+5);
         elseif  label(i+6)>0
            labelpred(i)=label(i+6);
         elseif  label(i+7)>0
            labelpred(i)=label(i+7);
         elseif  label(i+8)>0
            labelpred(i)=label(i+8);
        end
        
        
        
        % if label itself is apnea, we directly delete it 
        % . So use NAN
        if label(i)>0 
            labelpred(i)=nan;
        end
    end
end

end