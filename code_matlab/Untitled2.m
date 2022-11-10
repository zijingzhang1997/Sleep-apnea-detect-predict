label=zeros(length(TTime),1);
for i=1:size(annotations,1)
    for j=1:10
    if contains(stage(i),stg(j))
        ind=[round(RowTimes(i)):round(RowTimes(i))+round(duration(i))];
        label(ind+1)=lb(j);
        if label(ind+1)< lb(j)
             label(ind+1)=lb(j);
             
        end
        
    end
    end
end

[C,ia,ic] = unique(stage);
a_counts = accumarray(ic,1);
value_counts2 = [C, a_counts];


[C,ia,ic] = unique(label);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];