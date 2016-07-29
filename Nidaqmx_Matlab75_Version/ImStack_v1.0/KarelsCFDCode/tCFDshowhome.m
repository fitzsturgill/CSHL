
dir='d:\programs\images\';
base='h012a';
im=CFDshow(dir, base, 15,1);
%whos im
'select rectangle'
[x,y]=ginput(2)
x=sort(x);
y=sort(y);
rectangle('Position',[x(1),y(1),x(2)-x(1),y(2)-y(1)], 'EdgeColor','w')
'here1'
line=moviescan(im,x,y);
plot(line)
