function varargout = overlayGUI(varargin)
% OVERLAYGUI Application M-file for overlayGUI.fig
%    FIG = OVERLAYGUI launch overlayGUI GUI.
%    OVERLAYGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 03-Jan-2002 10:01:26

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

% --------------------------------------------------------------------
function varargout = generic_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.red1.
global gh state
genericCallback(h);


% --------------------------------------------------------------------
function varargout = red1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.red1.
global gh state
genericCallback(h);

set(gh.overlayGUI.red2, 'Value', 0);
genericCallback(gh.overlayGUI.red2);
set(gh.overlayGUI.red3, 'Value', 0);
genericCallback(gh.overlayGUI.red3);
set(gh.overlayGUI.green1, 'Value', 0);
genericCallback(gh.overlayGUI.green1);
set(gh.overlayGUI.blue1, 'Value', 0);
genericCallback(gh.overlayGUI.blue1);


% --------------------------------------------------------------------
function varargout = red2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.red2.
global gh state
genericCallback(h);

set(gh.overlayGUI.red1, 'Value', 0);
genericCallback(gh.overlayGUI.red1);
set(gh.overlayGUI.red3, 'Value', 0);
genericCallback(gh.overlayGUI.red3);
set(gh.overlayGUI.green2, 'Value', 0);
genericCallback(gh.overlayGUI.green2);
set(gh.overlayGUI.blue2, 'Value', 0);
genericCallback(gh.overlayGUI.blue2);



% --------------------------------------------------------------------
function varargout = red3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.red3.
global gh state
genericCallback(h);

set(gh.overlayGUI.red2, 'Value', 0);
genericCallback(gh.overlayGUI.red2);
set(gh.overlayGUI.red1, 'Value', 0);
genericCallback(gh.overlayGUI.red1);
set(gh.overlayGUI.green3, 'Value', 0);
genericCallback(gh.overlayGUI.green3);
set(gh.overlayGUI.blue3, 'Value', 0);
genericCallback(gh.overlayGUI.blue3);



% --------------------------------------------------------------------
function varargout = green1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.green1.
global gh state
genericCallback(h);

set(gh.overlayGUI.red1, 'Value', 0);
genericCallback(gh.overlayGUI.red1);
set(gh.overlayGUI.green3, 'Value', 0);
genericCallback(gh.overlayGUI.green3);
set(gh.overlayGUI.green2, 'Value', 0);
genericCallback(gh.overlayGUI.green2);
set(gh.overlayGUI.blue1, 'Value', 0);
genericCallback(gh.overlayGUI.blue1);



% --------------------------------------------------------------------
function varargout = green2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.green2.
global gh state
genericCallback(h);

set(gh.overlayGUI.green1, 'Value', 0);
genericCallback(gh.overlayGUI.green1);
set(gh.overlayGUI.green3, 'Value', 0);
genericCallback(gh.overlayGUI.green3);
set(gh.overlayGUI.red2, 'Value', 0);
genericCallback(gh.overlayGUI.red2);
set(gh.overlayGUI.blue2, 'Value', 0);
genericCallback(gh.overlayGUI.blue2);


% --------------------------------------------------------------------
function varargout = green3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.green3.
global gh state
genericCallback(h);

set(gh.overlayGUI.green2, 'Value', 0);
genericCallback(gh.overlayGUI.green2);
set(gh.overlayGUI.green1, 'Value', 0);
genericCallback(gh.overlayGUI.green1);
set(gh.overlayGUI.red3, 'Value', 0);
genericCallback(gh.overlayGUI.red3);
set(gh.overlayGUI.blue3, 'Value', 0);
genericCallback(gh.overlayGUI.blue3);


% --------------------------------------------------------------------
function varargout = blue1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.blue1.
global gh state
genericCallback(h);

set(gh.overlayGUI.red1, 'Value', 0);
genericCallback(gh.overlayGUI.red1);
set(gh.overlayGUI.blue3, 'Value', 0);
genericCallback(gh.overlayGUI.blue3);
set(gh.overlayGUI.blue2, 'Value', 0);
genericCallback(gh.overlayGUI.blue2);
set(gh.overlayGUI.green1, 'Value', 0);
genericCallback(gh.overlayGUI.green1);


% --------------------------------------------------------------------
function varargout = blue2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.blue2.
global gh state
genericCallback(h);

set(gh.overlayGUI.blue1, 'Value', 0);
genericCallback(gh.overlayGUI.blue1);
set(gh.overlayGUI.blue3, 'Value', 0);
genericCallback(gh.overlayGUI.blue3);
set(gh.overlayGUI.green2, 'Value', 0);
genericCallback(gh.overlayGUI.green2);
set(gh.overlayGUI.red2, 'Value', 0);
genericCallback(gh.overlayGUI.red2);


% --------------------------------------------------------------------
function varargout = blue3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.blue3.
global gh state
genericCallback(h);

set(gh.overlayGUI.blue2, 'Value', 0);
genericCallback(gh.overlayGUI.blue2);
set(gh.overlayGUI.blue1, 'Value', 0);
genericCallback(gh.overlayGUI.blue1);
set(gh.overlayGUI.green3, 'Value', 0);
genericCallback(gh.overlayGUI.green3);
set(gh.overlayGUI.red3, 'Value', 0);
genericCallback(gh.overlayGUI.red3);




% --------------------------------------------------------------------
function varargout = fileName1_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = fileName2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = fileName3_Callback(h, eventdata, handles, varargin)

