function SEM = nanSEM(A)

 
    SEM = std(A, 'omitnan') ./ sqrt(sum(~isnan(A), 1));