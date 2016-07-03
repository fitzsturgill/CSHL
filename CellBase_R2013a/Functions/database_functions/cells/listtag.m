function  list = listtag(xstr)
%LISTTAG   List CellBase tags.
%   L = LISTTAG(STR) returns a cell array containing the specified tag from
%   CellBase. Options for STR:
%   'properties', 'analysis', 'cell', 'rat'/'animal', 'session', 'tetrode'.
%   'properties' and 'analysis' are searched in ANALYSES, while others are
%   searched in CELLIDLIST.
%
%   See also LISTDIR and LISTFILES.

%   Edit log: BH 3/21/11, 5/3/12

% load cellbase
load(getpref('cellbase','fname'));

% Get tagged list
xstr = lower(char(xstr));
list = {};
if strncmp(xstr,'prop',4)
    n = 1;
    for i = 1:length(ANALYSES)
        for j = 1:length(ANALYSES(i).propnames)
            list{n} = char(ANALYSES(i).propnames(j));
            n = n + 1;
        end
    end
elseif strncmp(xstr,'anal',4)
    for i = 1:length(ANALYSES)
        list{i} = func2str(ANALYSES(i).funhandle);
    end
elseif strncmp(xstr,'rat',3) || strncmp(xstr,'animal',3)
    clist = char(CELLIDLIST);
    for i = 1:length(clist)
        [rat,remain] = strtok(clist(i,:),'_');
        list{i} = rat;
    end
    list = unique(list);
elseif strncmp(xstr,'ses',3)
    clist = char(CELLIDLIST);
    for i = 1:length(clist)
        [rat,remain] = strtok(clist(i,:),'_');
        [session,remain] = strtok(remain(2:end),'_');
        list{i,1} = rat(:)';
        list{i,2} = session(:)';
    end
    list = unique_cell(list);
elseif strncmp(xstr,'tetrode',3)
    clist = char(CELLIDLIST);
    for i = 1:length(clist)
        [rat,remain] = strtok(clist(i,:),'_');
        [session,remain] = strtok(remain(2:end),'_');
        [tetrode,remain] = strtok(remain(2:end),'.');
        list{i,1} = rat;
        list{i,2} = session;
        list{i,3} = tetrode;
    end
    list = unique_cell(list);
elseif strncmp(xstr,'cell',3)
    list = CELLIDLIST;
else
    disp('LISTTAG: unrecognized option.')
end