function varargout = standardModeGUI(varargin)
% STANDARDMODEGUI Application M-file for standardModeGUI.fig
%    FIG = STANDARDMODEGUI launch simpleModeGUI GUI.
%    STANDARDMODEGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 15-Jan-2001 15:14:20

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.

% --------------------------------------------------------------------
function varargout = generic_Callback(h, eventdata, handles, varargin)
% Stub for Callback of most uicontrol handles
	genericCallback(h);
	global state
	state.internal.secondsCounter=state.cycle.nextTimeDelay;
	updateGuiByGlobal('state.internal.secondsCounter');

function varargout = ZSlice_Callback(h, eventdata, handles, varargin)
% Stub for Callback of most uicontrol handles
	genericCallback(h);
	global state
	preallocateMemory;
	
function varargout = numberOfFrames_Callback(h, eventdata, handles, varargin)
% Stub for Callback of most uicontrol handles
	genericCallback(h);
	global state

	preAllocateMemory;
	
	if state.acq.dualLaserMode==1 % if the lasers are on simulataneously then nothing special
		sampleFactor=1;
	elseif state.acq.dualLaserMode==2
		sampleFactor=2;	% if they are alternating, then double the number of acqs before trigger the trigger function
	else
		disp('	setupInputChannels needs more for lasermodes');
	end
	set(state.daq.grabInput, 'SamplesPerTrigger', sampleFactor*state.internal.samplesPerFrame*state.acq.numberOfFrames);

	state.internal.needNewPcellRepeatedOutput=1;
	state.internal.needNewRepeatedMirrorOutput=1;
	applyChangesToOutput;
	
	% GRAB output: set number of frames in GRAB output object to drive mirrors


function varargout = averaging_Callback(h, eventdata, handles, varargin)
% Stub for Callback of most uicontrol handles
	genericCallback(h);
	global state
	preAllocateMemory;
	
% --------------------------------------------------------------------
function varargout = focusButton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.focusButton.
 	global state gh
	figure(gh.standardModeGUI.figure1);
	executeFocusCallback;	
		
% --------------------------------------------------------------------
function varargout = grabOneButton_Callback(h, eventdata, handles, varargin)
 	global state gh
	figure(gh.standardModeGUI.figure1);
	executeGrabOneCallback;

function varargout = dualLaserMode_Callback(h, eventdata, handles, varargin)
	genericCallback(h);
	updateDataForConfiguration;
