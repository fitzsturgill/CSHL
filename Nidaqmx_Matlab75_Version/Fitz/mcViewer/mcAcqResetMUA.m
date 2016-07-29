
function mcAcqResetMUA
    global state mcData
    
    
    
    includeChannels = find(mcChannelFieldToVector('MUAInclude'));
    if isempty(includeChannels)
        return
    end
    
    
    nCycles = length(state.cycle.delayList);
    
    
    for i = 1:nCycles
        wnAll = MUAAvgName('All', i);
        resetAverage(wnAll)
        for j = 1:length(includeChannels)
            channel = includeChannels(j);
            wn = MUAAvgName(channel, i);
            resetAverage(wn);
        end
    end
    