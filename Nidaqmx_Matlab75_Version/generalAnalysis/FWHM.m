function out=FWHM(data, offset)
	if nargin==1
		offset=min(data);
	end
	waveo('FWHMdata', data);
	[pd, py]=findpeaks(data, 1, offset, 0.3);
	waveo('FWHMx', pd-1);
	waveo('FWHMy', py);
	out=pd(2)-pd(1);
	