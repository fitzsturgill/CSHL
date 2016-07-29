function [h,rectPos]=drawSolidRectangle
% Thsi function creates a rubber bix
% return the handle to the rectangle and the position of the rectangle.

pos=getrect(gca);
rectPos=round(pos);
h=rectangle('Position',round(pos),'edgecolor',[1 0 0]);