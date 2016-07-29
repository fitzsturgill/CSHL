function out = showSignificantPairsN(c, type)
global gh state
out = [];
% Takes output from multcompare and displays significant pairs.

if nargin < 2
	type = '';
end

totalInteractions = size(c,1);
counter = 1;
for i = 1:totalInteractions
if (c(i,3) < 0 & c(i,5) < 0) | (c(i,3) >= 0 & c(i,5) > 0) 
	disp(['Pair ' num2str(c(i,1)) ' and ' num2str(c(i,2)) ' are significantly different in ' type '.']);
	out(counter, 1:2) = c(i,1:2);
	counter = counter+1;
end
end
if isempty(out)
	disp(['No significant pairs found for ' type]);
end



		