function markBaseline(number)    

    global state
	global gh
	
    if nargin<1
        number=state.notebook.notebookNumber;
    end
    
    current=clock;
% 	if state.epoch>1
		state.epoch=state.epoch+1;
		updateGUIbyGlobal('state.epoch');
% 	end
	genericCallback(gh.timerMainControls.epoch);
	addEntryToNotebook(number, ...
        ['baseline' ':' '000' ':' num2str(etime(current, state.internal.startupTime)) ':epoch:' num2str(state.epoch)]);
	disp(['*** ' 'Baseline ' num2str(etime(current, state.internal.startupTime)) 'sec from startup***'])

    
    % 		[datestr(clock,13) ', Epoch ' num2str(state.epoch) ', Acq ' num2str(state.files.lastAcquisition) ' : ' state.notebook.newEntry]);
	closeNotebookEntry;

	timerLoop;