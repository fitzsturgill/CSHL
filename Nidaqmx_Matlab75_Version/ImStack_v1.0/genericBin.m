function out = genericBin(image, binX, binY, binZ);
global gh state

% function bins the image by binX, binY, and binZ...
type = class(image);
image = double(image);

Ny = size(image,1); 		% Tells function the number of rows of image
Nx = size(image,2);			% Tells function the number of columns of image

dimImage = ndims(image);		% Tells function the dimensions of image as a scalar
if dimImage == 3				
	Nz = size(image,3);				% Tells function the number of frames of image
end

if binX == 1 & binY == 1 & binZ == 1
	out = image;
	return
end

if dimImage > 2
	
	if binX == 1 & binY > 1 & binZ > 1
		image = permute(image, [1 3 2]);									% turn frames into columns
		image = reshape(image, binY, (Ny/binY), binZ, (Nz/binZ), Nx);		% Reshapes Array into 5 dimensions 
		image = permute(image, [1 3 2 4 5]);								% Interchanges columsn and frames
		image = reshape(image, (binY*binZ),((Ny*Nz*Nx)/(binY*binZ)));		% Puts binning along the rows and columns	
		image = sum(image);
		image = (reshape(image, Ny/binY, Nz/binZ, Nx))/(binY*binZ);
		eval(['out = ' type '(permute(image, [1 3 2]));']);	
	end
	
	if binY == 1 & binX > 1 & binZ > 1
		image = permute(image, [3 2 1]);									% turn frames into columns
		image = reshape(image, binZ, (Nz/binZ), binX, (Nx/binX), Ny);		% Reshapes Array into 5 dimensions 
		image = permute(image, [1 3 2 4 5]);								% Interchanges columsn and frames
		image = reshape(image, (binX*binZ),((Ny*Nz*Nx)/(binZ*binX)));		% Puts binning along the rows and columns	
		image = sum(image);
		image = (reshape(image, Nz/binZ, Nx/binX, Nx))/(binX*binZ);
		eval(['out = ' type '(permute(image, [3 2 1]));']);	
	end
	
	if binZ == 1 & binY > 1 & binX > 1
		image = reshape(image, binY, (Ny/binY), binX, (Nx/binX), Nz);		% Reshapes Array into 5 dimensions 
		image = permute(image, [1 3 2 4 5]);								% Interchanges columsn and frames
		image = reshape(image, (binY*binX),((Nx*Ny*Nz)/(binY*binX)));		% Puts binning along the rows and columns		
		image = sum(image);													% Sums accross all rows	
		%out = (reshape(image, (Ny/binY), (Nx/binX), Nz))/(binX*binY);		% Reshapes matrix into its new dimensions and normalizes
		eval(['out = ' type '((reshape(image, (Ny/binY), (Nx/binX), Nz))/(binX*binY));']);	
	end
	
	if binY > 1 & binX > 1 & binZ > 1
		image = reshape(image, binY, (Ny/binY), binX, (Nx/binX), binZ, (Nz/binZ));	% Reshapes Array into 6dimensions 
		image = permute(image, [1 3 5 2 4 6]);										% Interchanges columsn and frames
		image = reshape(image, (binY*binX*binZ),((Nx*Ny*Nz)/(binY*binX*binZ)));		% Puts binning along the rows and columns		
		image = sum(image);															% Sums accross all rows
	
		if ndims(image) == 2
			%out = (reshape(image, (Ny/binY), (Nx/binX)))/(binX*binY*binZ);					% Reshpapes into new form
			eval(['out = ' type '((reshape(image, (Ny/binY), (Nx/binX)))/(binX*binY*binZ));']);	
		else
			out = (reshape(image, (Ny/binY), (Nx/binX), (Nz/binZ)))/(binX*binY*binZ);			% Reshpapes into new form
			eval(['out = ' type '((reshape(image, (Ny/binY), (Nx/binX), (Nz/binZ)))/(binX*binY*binZ));']);	
		end
	end
	
	if binX > 1 & binY == 1 & binZ == 1
		image = permute(image, [2 1 3]);
		image = reshape(image, binX, Nx/binX, Ny,Nz);
		image = sum(image);
		image = reshape(image, Nx/binX,Ny,Nz)/(binX);
		%out = permute(image, [2 1 3]);
		eval(['out = ' type '(permute(image, [2 1 3]));']);	
	end
	
	if binY > 1 & binX == 1 & binZ == 1
		image = reshape(image, binY, Ny/binY, Nx,Nz);
		image = sum(image);
		image = reshape(image, Ny/binY,Nx,Nz)/(binY);
		eval(['out = ' type '(reshape(image, Ny/binY,Nx,Nz)/(binY));']);
	end
	
	if binZ > 1 & binY == 1 & binX == 1
		image = permute(image, [3 2 1]);
		image = reshape(image, (binZ), Nz/binZ, Nx,Ny);	% Reshapes Array into 6dimensions 
		image = sum(image);		
		image = squeeze(image);
		image = permute(image, [3 2 1])/(binZ);
		image = squeeze(image);
		% Sums accross all rows

		if ndims(image) == 2
			%out = (reshape(image, (Ny/binY), (Nx/binX)))/(binX*binY*binZ);					% Reshpapes into new form
			eval(['out = ' type '(image);']);	
		else
			%out = (reshape(image, (Ny/binY), (Nx/binX), (Nz/binZ)))/(binX*binY*binZ);			% Reshpapes into new form
			eval(['out = ' type '(image);']);	
		end
		
	end
	
else
	image = reshape(image, binY, (Ny/binY), binX, (Nx/binX));		% Reshapes Array into 5 dimensions 
	image = permute(image, [1 3 2 4]);								% Interchanges columsn and frames
	image = reshape(image, (binY*binX),((Nx*Ny)/(binY*binX)));		% Puts binning along the rows and columns		
	image = sum(image);													% Sums accross all rows	
	out = (reshape(image, (Ny/binY), (Nx/binX)))/(binX*binY);		% Reshapes matrix into its new dimensions and normalizes
	eval(['out = ' type '((reshape(image, (Ny/binY), (Nx/binX)))/(binX*binY));']);	
end

		

		
		
	