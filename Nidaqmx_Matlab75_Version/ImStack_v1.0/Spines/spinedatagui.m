function varargout = SpineDataGUI(varargin)
% SPINEDATAGUI Application M-file for SpineDataGUI.fig
%    FIG = SPINEDATAGUI launch SpineDataGUI GUI.
%    SPINEDATAGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 27-Aug-2001 14:14:56

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
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
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
function varargout = overallDensErr_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);



% --------------------------------------------------------------------
function varargout = overallDensity3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = overallDenErr3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = meanDensity3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = meanDenErr3_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = analyzeByFile_Callback(h, eventdata, handles, varargin)

global gh state
genericCallback(h);
state.imageProc.spineData.analyzeByNeuron=0;
updateGUIByGlobal('state.imageProc.spineData.analyzeByNeuron');

% --------------------------------------------------------------------
function varargout = analyzeByNeuron_Callback(h, eventdata, handles, varargin)

global gh state
genericCallback(h);
state.imageProc.spineData.analyzeByFile=0;
updateGUIByGlobal('state.imageProc.spineData.analyzeByFile');




% --------------------------------------------------------------------
function varargout = NumberOfNeurons_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);



% --------------------------------------------------------------------
function varargout = Untitled_16_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_19_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_20_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_21_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_28_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = add_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.add.
addSpineDataToList;


% --------------------------------------------------------------------
function varargout = remove_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.remove.
clearSpineDataStats;



% --------------------------------------------------------------------
function varargout = anovaP_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.anovaP.



% --------------------------------------------------------------------
function varargout = annovaP_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.annovaP.
genericCallback(h);


% --------------------------------------------------------------------
function varargout = Untitled_39_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = annova2d_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.annova2d.
genericCallback(h);


% --------------------------------------------------------------------
function varargout = annova3d_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.annova3d.
genericCallback(h);


% --------------------------------------------------------------------
function varargout = annovaLen_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.annovaLen.
genericCallback(h);


% --------------------------------------------------------------------
function varargout = annovaVol_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.annovaVol.
genericCallback(h);


% --------------------------------------------------------------------
function varargout = numberOfSpines_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_40_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_41_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_42_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_43_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_44_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_46_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_47_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = Untitled_48_Callback(h, eventdata, handles, varargin)

