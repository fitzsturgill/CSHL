function out=GetSweepList(directory, base, extension)
% GetSweepList searches Directory for files of type 'BaseSweep.Extension'
%		GetSweepList returns an array of Sweep numbers 
%		that satisfy the convention 'Base###.Extension'
% 		Note: extensions are of type '.ext'
%		
%   Karel Svoboda 1/16/00 Matlab 5.3
%	 svoboda@cshl.org
%
%
if nargin < 3 
     disp(['-------------------------------']);  
     disp([' GetSweepList v',version])
     disp(['-------------------------------']);
     disp([' usage: GetSweepList(directory, base, extension)']);
     error(['### not enough parameters to proceed']); 
  end
 
nam=strcat(directory, base,'*',extension);
flist=dir(nam);

len=length(flist);
numlist=1:len;

for i=1:len
   cc=strtok(flist(i).name,'.cfd');
   lencc=length(cc);
   cc=cc(lencc-3+1:lencc);
   numlist(i)=str2num(cc);
end

out=numlist;
