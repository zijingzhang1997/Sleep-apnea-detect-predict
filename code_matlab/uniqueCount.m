function label_counts=uniqueCount(x)

[C,ia,ic] = unique(x);
a_counts = accumarray(ic,1);
label_counts = [C, a_counts];
end