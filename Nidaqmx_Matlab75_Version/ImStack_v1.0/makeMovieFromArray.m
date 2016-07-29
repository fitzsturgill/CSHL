function makeMovieFromArray(array,repeat,fps)
global state gh 

frames = size(array,3);



% Record the movie

for j = 1:frames 
   	image(array(:,:,j));
	axis image;
	colormap(gray);
	F(j) = getframe;
end
