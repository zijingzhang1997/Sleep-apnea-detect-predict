% consider subject variation 
% substract baseline of each subject 
%baseline = mean of features in normal epoch

function [labelNew,featureNew]=feature_calibrate(label,feature,opt)



labelNew=label;
featureNew=feature;


if opt==1
    [idxNan,~]=find(feature~=feature);
    featureNew(idxNan,:)=[];
    labelNew(idxNan)=[];
    
    
  
    
    % baseline set as the mean of all normal epochs
    [idx,~]=find(labelNew==0);
    baseline=mean(featureNew(idx,:),1);
    
    
    
    for i=1:length(featureNew)
        featureNew(i,:)=featureNew(i,:)-baseline;
    end
    
end


end