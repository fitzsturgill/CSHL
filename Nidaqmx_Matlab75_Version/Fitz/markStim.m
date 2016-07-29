function markStim(length, number)    

    global state
	global gh
	
    if nargin<2
        number=state.notebook.notebookNumber;
    end
	
	if nargin<1
		length=60;
	end
    
    current=clock;
	
	state.epoch=state.epoch+1;
	updateGUIByGlobal('state.epoch');
	
	addEntryToNotebook(number, ...
        ['stim' ':' num2str(length) ':' num2str(etime(current, state.internal.startupTime)) ':epoch:' num2str(state.epoch)]);
	disp(['*** ' num2str(length) 'sec Puff at ' num2str(etime(current, state.internal.startupTime)) 'sec from startup***'])

    
    % 		[datestr(clock,13) ', Epoch ' num2str(state.epoch) ', Acq ' num2str(state.files.lastAcquisition) ' : ' state.notebook.newEntry]);

	closeNotebookEntry;

	
	genericCallback(gh.timerMainControls.epoch);
	timerLoop;
	
	state.epoch=state.epoch+1;
	updateGUIByGlobal('state.epoch');
