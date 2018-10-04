%     Script to delete weird spurious files that Bonsai spit out, there are
%     extra files for certain triggers that I can detect because their
%     datemodified timestamp is < 2 seconds from the previous timestamp
%     which is impossible given that Bpod trials are at least 11 seconds
%     long + the ITI

   
%     cd('Z:\FitzRig2\Data\DC_47\lickNoLick_Odor_v2\Session Data\Whisk_180621');
%     fs_whisk = dir('*hisk_*.avi');
%     fileList_whisk = {};
%     [fileList_whisk{1:length(fs_whisk)}] = fs_whisk(:).name; % see deal documentation I think...
%     [fileList_whisk,ix_whisk] = sort_nat(fileList_whisk); % alphanumeric sorting
%     dmDelta_whisk = seconds(diff(datetime({fs_whisk(ix_whisk).date})));
%     
%     todelete = fileList_whisk(find(dmDelta_whisk < 3) + 1);
%     for counter = 1:length(todelete)
%         eval(sprintf('delete %s', todelete{counter}));
%     end
%     
    cd('Z:\FitzRig2\Data\DC_46\lickNoLick_Odor_v2\Session Data\WhiskDiff_180614');
    fs_whiskDiff = dir('*hiskDiff_*.csv');
    fileList_whiskDiff = {};
    [fileList_whiskDiff{1:length(fs_whiskDiff)}] = fs_whiskDiff(:).name; % see deal documentation I think...
    [fileList_whiskDiff,ix_whiskDiff] = sort_nat(fileList_whiskDiff); % alphanumeric sorting
    dmDelta_whiskDiff = seconds(diff(datetime({fs_whiskDiff(ix_whiskDiff).date})));
    
    todelete = fileList_whiskDiff(find(dmDelta_whiskDiff < 3) + 1);
    for counter = 1:length(todelete)
        eval(sprintf('delete %s', todelete{counter}));
    end

    
    
    