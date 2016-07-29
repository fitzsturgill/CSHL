load mri
h = imshow(D(:,:,:,1), [1 100]); 
%set(h,'erasemode','xor'); 
        for i = 2:10 
           set(h,'cdata',D(:,:,:,i))
           pause
        end