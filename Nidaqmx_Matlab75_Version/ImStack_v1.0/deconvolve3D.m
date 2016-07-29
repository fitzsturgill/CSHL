function out=deconvolve3D(image,PSF,iterations)
% thsi function will take a 3D dimensional array and deconvolve it with the PSF
% IT will do this iterations times.
% This uses a standard wiener deconvolution algorithm.

out=[];

%image = edgetaper(image,PSF);   % Reduces end ringing
if nargin < 2 | nargin > 3
    error('deconvolve3D:Need to supply a PSF and an image.  MAx three inputs.');
elseif nargin == 2
    iterations=1;
end
frames=size(image,3);
h = waitbar(0,'Deconvolution of Multi-Tif', 'Name', 'Deconvolution 3D');
for i = 1:iterations
	for j=1:frames
		waitbar(j/frames, h, ['Iteration ' num2str(i) ' of ' num2str(iterations) '; Frame ' num2str(j) ' of ' num2str(frames)]);
    	image(:,:,j) = deconvwnr(image(:,:,j),double(PSF));
	end
end
close(h);
out=image;
