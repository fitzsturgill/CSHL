function varargout = fileCounterGUI(varargin)
% FILECOUNTERGUI Application M-file for fileCounterGUI.fig
%    FIG = FILECOUNTERGUI launch fileCounterGUI GUI.
%    FILECOUNTERGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 11-Apr-2001 10:45:02

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

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

function varargout = generic_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.
global state gh
value = get(gh.imageProcessingGUI.fileName, 'Value');
genericCallbackCell(h);

function varargout = genericNormal_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.
global state gh
genericCallback(h);


function varargout = loadImage_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.
global state gh

pathName = get(gh.fileCounterGUI.pathName, 'String');
baseName = get(gh.fileCounterGUI.baseName, 'String');

	
value = get(gh.fileCounterGUI.filetype, 'Value');
ext = get(gh.fileCounterGUI.filetype, 'String');
ext = ext{value};
acquisitionNumber = get(gh.fileCounterGUI.acquisitionCounter, 'String');

	if state.imageProc.appendZeros == 1
		if str2num(acquisitionNumber) < 10
			baseName = [baseName '0'];
		else
		end
	end

filename = [pathName baseName num2str(acquisitionNumber) ext];
	if strcmp(ext, '.tif')
		loadImageFromName(filename,pathName);
	elseif strcmp(ext, '.jpg')
		state.imageProc.currentjpeg = imread(filename);
		state.imageProc.currentjpeg = rgb2gray(state.imageProc.currentjpeg);
		loadImageFromArray('state.imageProc.currentjpeg');
    elseif strcmp(ext, '.CFD')
		seeGUI('gh.cfdGUI.figure1');
		set(gh.fileCounterGUI.filetype, 'Value', 3);
        openACFD(filename);
		hideGUI('gh.cfdGUI.figure1');
    elseif strcmp(ext, '.cfd')
		set(gh.fileCounterGUI.filetype, 'Value', 3);
        seeGUI('gh.cfdGUI.figure1');
        openACFD(filename);
		hideGUI('gh.cfdGUI.figure1');
	end
state.imageProc.acquisitionNumber = state.imageProc.acquisitionNumber+1;
updateGUIByGlobal('state.imageProc.acquisitionNumber');


% --------------------------------------------------------------------
function varargout = filetype_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.filetype.
