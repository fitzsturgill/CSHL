    function out=MUAAvgName(channel, cyclePos, epoch)
        global state
        
        if nargin < 1
            channel = 'ALL';
        elseif isnumeric(channel)
            channel = num2str(channel);
        end
        
        if nargin < 2
            cyclePos = num2str(state.cycle.currentCyclePosition);
        elseif isnumeric(cyclePos)
            cyclePos = num2str(cyclePos);
        end
        
        if nargin < 3
            epoch = num2str(state.epoch);
        elseif isnumeric(epoch)
            epoch= num2str(epoch);
        end


        



        out=[state.files.baseName 'e' epoch 'ch' channel '_cyc' cyclePos '_MUA'];