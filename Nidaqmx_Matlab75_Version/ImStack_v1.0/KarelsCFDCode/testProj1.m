
global head_info

fn='D:\data\brian\chronic_data\c083b001.cfd';
head_info=cf4header(fn);

for i=1:2;
   im0=ProjStack1(fn, 1, [1 25], 'x', 10*(i-1));
end   
imshow(im{1});
pause;
imshow(im{1});
