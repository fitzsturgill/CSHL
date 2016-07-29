function mcAcqInitializePhysDaq
    global state
    % FS MOD-   multichannel acquisition, add channels here as they won't
    % change
	if ~isempty(state.phys.mcAcq.mcInputBoardIndex) % should I have a variable mcAcqOn?
		state.phys.mcAcq.mcInputDevice = analoginput('nidaq', state.phys.mcAcq.mcInputBoardIndex);
		set(state.phys.mcAcq.mcInputDevice, 'TriggerType', 'HwDigital');		
        set(state.phys.mcAcq.mcInputDevice, 'InputType', 'SingleEnded');
        set(state.phys.mcAcq.mcInputDevice, 'ExternalTriggerDriveLine', 'RTSI0'); % generate a pulse on the RTSI0 channel when this engine is triggered...        
% 		set(state.phys.mcAcq.mcInputDevice, 'SampleRate', state.phys.mcAcq.mcInputRate * 1000);
        % use setverify in case the device can't be set to the specified
        % value
		state.phys.mcAcq.mcInputRate = setverify(state.phys.mcAcq.mcInputDevice, 'SampleRate', state.phys.mcAcq.mcInputRate * 1000) / 1000;
        updateGUIByGlobal('state.phys.mcAcq.mcInputRate');
%         set(state.phys.mcAcq.mcInputDevice, 'SamplesAcquiredFcn', {'processMCData'});
        disp('mcAcqInitializePhysDaq: mcInputDevice SamplesAcquiredFcn not set, it is called from processPhysData');

        
        % also create and fill channel structure to control how acquired
        % data is displayed
        channelStruct = struct(...
            'Name', '', ...   
            'Show', 0, ...
            'ShowFilter', 0, ...
            'MUAInclude', 0, ...
            'LowPass', 0, ...
            'HighPass', 0 ...
        );
        % initialize channel structure
        state.phys.mcAcq.channel = repmat(channelStruct, 1, state.phys.mcAcq.totalChannels); % create channel structure to be modified                    
        for i = 1:state.phys.mcAcq.totalChannels
            % add channels
            channel = addchannel(state.phys.mcAcq.mcInputDevice, i-1, i); % hwchannel is 0 based, right???? double check...
            channel.InputRange = [-10 10];
            channel.UnitsRange = [-10 10];
%             set(state.phys.mcAcq.mcChannels(i), 'InputRange', [-10 10]);
%             set(state.phys.mcAcq.mcChannels(i), 'UnitsRange', [-10 10]);

            % further initialization of channel structure
            % use channelNames to provide channel names if possible,
            % otherwise set default names and define channelNames
            if i<=state.phys.mcAcq.mcNChannels
%                 state.phys.mcAcq.channel(i).Name = ['mc' num2str(i)];
%                 state.phys.mcAcq.mcChannelNames{1, i} = ['mc' num2str(i)];
                state.phys.mcAcq.channel(i).Show = 1;
            end           
        end
        mcAcqUpdateChannelNames; % use channel names based upon channel ordering
    state.phys.mcAcq.displayData=zeros(1, state.phys.mcAcq.totalChannels);
    state.phys.mcAcq.displayXData=0;
    state.phys.mcAcq.displayThreshData=zeros(1, state.phys.mcAcq.totalChannels);
    end
    % end MOD