function y = decimate2(x, r, dim)
%     Accepts 2d matrices, treats columns as independent vectors
% dim- optional, specify to decimate along columns or rows

    if nargin < 3
        dim = 1;
    end
   
    if dim == 2
        x = x';
    end
    y = zeros(ceil(size(x,1)/r), size(x, 2));
    if isvector(x)
        y = decimate(x, r);
        return
    else
        for i=1:size(x, 2)
            y(:,i) = decimate(x(:,i), r);
        end
    end
    
    if dim == 2
        y = y'; % switch it back
    end