function out = morphAnovaN(data,group1, group2,ss,intxns,names)

totalGroup= {group1' group2'};
[p2d,tbl2d,stat2d] = anovan(data',totalGroup,intxns,ss,names);
out = p2d;