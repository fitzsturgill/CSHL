function mcAcqOlf_callback(obj, event)
    % this function is triggered by timer property of DIO object olfDEvice
    % it switches on the odor by actuating the valves.  
    % it then changes the period of the timer and restarts the timer in
    % order to switch after the odor after the olfDelay
    global state
    
    if ~state.phys.mcAcq.olfShuntEnabled
        if ~state.phys.mcAcq.olfValveStatus % odor is off
            mcAcqValveSwitch(state.cycle.mcOlfValve);  % turn valves on
            state.phys.mcAcq.olfValveStatus = 1;
            stop(state.phys.mcAcq.olfDevice);
            set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.cycle.mcOlfDuration); % reset timer period
            start(state.phys.mcAcq.olfDevice);
            disp('odor turned on');
%             %%
%             Fs = 1000;      %# Samples per second
%             toneFreq = 50;  %# Tone frequency, in Hertz
%             nSeconds = 2;   %# Duration of the sound
%             y = sin(linspace(0, nSeconds*toneFreq*2*pi, round(nSeconds*Fs)));
% 
%             When played at 1 kHz using the SOUND function, this vector will generate a 50 Hz tone for 2 seconds:
% 
%             sound(y, Fs);  %# Play sound at sampling rate Fs            
            %%
        else % odor is on
            mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve);
            stop(state.phys.mcAcq.olfDevice);
            state.phys.mcAcq.olfValveStatus = 0;
            disp('odor turned off');
        end
    else
        if state.phys.mcAcq.olfValveStatus == 0 % odor is off
            mcAcqValveSwitch(state.cycle.mcOlfValve, 1);  % turn valves on, begin shunt
            state.phys.mcAcq.olfValveStatus = - 1;  % indicates shunting in progress
            stop(state.phys.mcAcq.olfDevice);
            set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.phys.mcAcq.olfShuntDuration); % reset timer period
            start(state.phys.mcAcq.olfDevice);
            disp('odor turned on, shunting');
        elseif state.phys.mcAcq.olfValveStatus == -1 % shunting is in progress
            mcAcqValveSwitch(state.cycle.mcOlfValve, 0);  % leave valves on, end shunt
%             mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve, 1);  % leave valves on, end shunt
            state.phys.mcAcq.olfValveStatus = 1;  % indicates odor is on
            stop(state.phys.mcAcq.olfDevice);
            set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.cycle.mcOlfDuration); % reset timer period
            start(state.phys.mcAcq.olfDevice);
            disp('odor turned on');
        else % odor is on
            mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve, 0);
            stop(state.phys.mcAcq.olfDevice);
            state.phys.mcAcq.olfValveStatus = 0;
            disp('odor turned off');
        end
    end
    
    
    
    
%     if ~state.phys.mcAcq.olfValveStatus % odor is off
%         mcAcqValveSwitch(state.cycle.mcOlfValve);  % turn valves on
%         state.phys.mcAcq.olfValveStatus = 1;
%         stop(state.phys.mcAcq.olfDevice);
%         set(state.phys.mcAcq.olfDevice, 'TimerPeriod', state.cycle.mcOlfDuration); % reset timer period
%         start(state.phys.mcAcq.olfDevice);
%         disp('odor turned on');
%     else % odor is on
%         mcAcqValveSwitch(state.phys.mcAcq.olfDefaultValve);
%         stop(state.phys.mcAcq.olfDevice);
%         state.phys.mcAcq.olfValveStatus = 0;
%         disp('odor turned off');
%     end    