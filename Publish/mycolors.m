function rgb = mycolors(condition)
    if nargin < 1
        condition = '';
    end
    color_table = {...
        'chat', [171 55 214]/256;... % purple
        'dat', [237 125 49]/256;...
        'licks', [0 0 0];...
        'shock', [128/255 255/255 0];... 
        };
   
    row = find(strcmpi(color_table(:,1), condition));
    if row
        rgb = color_table{row, 2};
    else
        disp('choices:');
        disp(color_table(:,1));
    end
    