function varargout = mathGUI(varargin)
% MATHGUI Application M-file for mathGUI.fig
%    FIG = MATHGUI launch mathGUI GUI.
%    MATHGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 15-Feb-2002 15:34:05

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
function varargout = fileName1_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = startFile1_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = startFile1Slider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = fileName2_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = startFile2_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = endFile1_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = endFile1Slider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = endFile2_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = endFile2Slider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = operation_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);
str=get(h,'String');
val=get(h,'Value');
str=str{val};
if strcmp(str,'Register')
	set([gh.mathGUI.register gh.mathGUI.lut gh.mathGUI.luttext gh.mathGUI.transform gh.mathGUI.transformtext],'visible','on');
else
	set([ gh.mathGUI.lut gh.mathGUI.luttext gh.mathGUI.register gh.mathGUI.transform gh.mathGUI.transformtext],'visible','off');
end

% --------------------------------------------------------------------
function varargout = weight_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = loadImage_Callback(h, eventdata, handles, varargin)
loadMathImage;




% --------------------------------------------------------------------
function varargout = startFile2Slider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = register_Callback(h, eventdata, handles, varargin)
global gh state input_points base_points mytform

% Do the registering and transform....

set(gh.mathGUI.figure1,'Pointer', 'watch');
val=get(gh.mathGUI.transform, 'Value');
str=get(gh.mathGUI.transform, 'String');
transform=str{val};
input_pts_adj= cpcorr(input_points, base_points, state.imageProc.mathGUI.inputImage, state.imageProc.mathGUI.baseImage);
mytform = cp2tform(input_pts_adj,base_points,transform);

[rows cols]=size(state.imageProc.mathGUI.inputImage);
try
	state.imageProc.mathGUI.modImage = imtransform(state.imageProc.mathGUI.inputImage,mytform, 'Xdata',[1 cols], 'Ydata',[1 rows]);
	loadImageFromArray('state.imageProc.mathGUI.modImage');
	set(gh.mathGUI.figure1,'Pointer', 'Arrow');
catch
	beep;
	disp('Unable to register image. Select Mopre poitns or try a different transform');
	set(gh.mathGUI.figure1,'Pointer', 'Arrow');
end



% --------------------------------------------------------------------
function varargout = add_Callback(h, eventdata, handles, varargin)
addCurrentImgToStack;



% --------------------------------------------------------------------
function varargout = remove_Callback(h, eventdata, handles, varargin)
global state
q=size(state.imageProc.imgStack,3);
if q == 1
	state.imageProc.imgStack = [];
else
	state.imageProc.imgStack =state.imageProc.imgStack(:,:,1:end-1);
end


% --------------------------------------------------------------------
function varargout = clear_Callback(h, eventdata, handles, varargin)
global state

state.imageProc.imgStack = [];



% --------------------------------------------------------------------
function varargout = load_Callback(h, eventdata, handles, varargin)
global state

if ~isempty(state.imageProc.imgStack)
	loadImageFromArray('state.imageProc.imgStack');
else
	disp('Cant load an image that is empty');
end


% --------------------------------------------------------------------
function varargout = transform_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = lut_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = mathValue_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = useValue_Callback(h, eventdata, handles, varargin)
genericCallback(h);
