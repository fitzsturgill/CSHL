function lengths = computeLengthOfLine(line)
global state

total = length(line);
for i = 1:total
	Xcoords = get(line(i), 'XData');
	Ycoords = get(line(i), 'YData');
	Xcoordsize= length(Xcoords);
	part = [];
	for j = 1:Xcoordsize-1
		part(j) = ((state.imageProc.spine.micronsperpixelX*(Xcoords(j)-Xcoords(j+1)))^2 + (state.imageProc.spine.micronsperpixelY*(Ycoords(j)-Ycoords(j+1)))^2)^.5;
	end
		lengths(i) = sum(part);
end
