function out=NextFilename(target_directory, filename, ext)
% CFDread2 generaes an unoccupied filename 
% this function depends on the naming convention a001a001.ext.
%		
%   Karel Svoboda 1/10/00 Matlab 5.3
%	 svoboda@cshl.org
%


len=length(filename);

cc=filename((len-11):(len-4));
ccc=strcat(target_directory,cc,ext);
fid = fopen(ccc,'r');
i=1;
while fid > 0
   num=num2str(i)
   ccc=strcat(target_directory,cc,'_',num,ext);
   fid = fopen(ccc,'r');
end

out=ccc;