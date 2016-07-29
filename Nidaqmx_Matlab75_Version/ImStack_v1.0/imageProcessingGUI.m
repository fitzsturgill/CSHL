function varargout = imageProcessingGUI(varargin)
% IMAGEPROCESSINGGUI Application M-file for imageProcessingGUI.fig
%    FIG = IMAGEPROCESSINGGUI launch imageProcessingGUI GUI.
%    IMAGEPROCESSINGGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 07-Jun-2002 13:57:01

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

function varargout = fileName_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrame.
global gh state
genericCallback(h);
value = get(h, 'Value');
switchFileName(value);


% --------------------------------------------------------------------
function varargout = lowpixel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrame.
global gh state
genericCallbackCell(h);
value = get(gh.imageProcessingGUI.fileName, 'Value');
set(state.imageProc.internal.axis{value}, 'CLim', [state.imageProc.lowPixelValue ...
		state.imageProc.highPixelValue]);
if state.imageProc.updateLUTMax 
	
	set([gca gh.maxProjectionGUI.xyaxis gh.maxProjectionGUI.xzaxis gh.maxProjectionGUI.yzaxis ], 'CLim', [state.imageProc.lowPixelValue ...
			state.imageProc.highPixelValue]);
end

% --------------------------------------------------------------------
function varargout = highpixel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrame.
global gh state
genericCallbackCell(h);
value = get(gh.imageProcessingGUI.fileName, 'Value');
set(state.imageProc.internal.axis{value}, 'CLim', [state.imageProc.lowPixelValue ...
		state.imageProc.highPixelValue]);

if state.imageProc.updateLUTMax 
	set([gca gh.maxProjectionGUI.xyaxis gh.maxProjectionGUI.xzaxis gh.maxProjectionGUI.yzaxis ], 'CLim', [state.imageProc.lowPixelValue ...
			state.imageProc.highPixelValue]);
end



% --------------------------------------------------------------------
function varargout = currentFrame_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrame.
global state gh

value = get(gh.imageProcessingGUI.fileName, 'Value');
genericCallbackCell(h);
showCurrentFrame(value);
state.imageProc.internal.meanIntensity = mean2(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame));
updateGUIByGlobal('state.imageProc.internal.meanIntensity');
state.imageProc.internal.sumIntensity = sum(sum(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame)));
updateGUIByGlobal('state.imageProc.internal.sumIntensity');
if state.imageProc.autoIntensity == 1
	state.imageProc.lowPixelValue = min(min(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame)));
	state.imageProc.highPixelValue = (.01*state.imageProc.LUTpercentmax)*double(max(max(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame))));
	updateGUIByGlobal('state.imageProc.lowPixelValue');
	updateGUIByGlobal('state.imageProc.highPixelValue');
	set(state.imageProc.internal.axis{value}, 'CLim', [state.imageProc.lowPixelValue ...
		state.imageProc.highPixelValue]);
end
	


% --------------------------------------------------------------------
function varargout = currentFrameSlider_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
value = get(gh.imageProcessingGUI.fileName, 'Value');
genericCallbackCell(h);
showCurrentFrame(value);
state.imageProc.internal.meanIntensity = mean2(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame));
updateGUIByGlobal('state.imageProc.internal.meanIntensity');
state.imageProc.internal.sumIntensity = sum(sum(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame)));
updateGUIByGlobal('state.imageProc.internal.sumIntensity');
if state.imageProc.autoIntensity == 1
	state.imageProc.lowPixelValue = min(min(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame)));
	state.imageProc.highPixelValue = (.01*state.imageProc.LUTpercentmax)*double(max(max(state.imageProc.cell.currentImage{value}(:,:,state.imageProc.currentFrame))));
	updateGUIByGlobal('state.imageProc.lowPixelValue');
	updateGUIByGlobal('state.imageProc.highPixelValue');
	set(state.imageProc.internal.axis{value}, 'CLim', [state.imageProc.lowPixelValue ...
		state.imageProc.highPixelValue]);
end



% --------------------------------------------------------------------
function varargout = totalFramesSlider_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.totalFrames.
global state gh
genericCallbackCell(h);

% --------------------------------------------------------------------
function varargout = generic_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallbackCell(h);

% --------------------------------------------------------------------
function varargout = updateClim_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.updateClim.
global gh state
genericCallback(h);


% --------------------------------------------------------------------
function varargout = meanIntensity_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.meanIntensity.
global state gh





% --------------------------------------------------------------------
function varargout = Untitled_1_Callback(h, eventdata, handles, varargin)
global state gh
genericCallback(h);



% --------------------------------------------------------------------
function varargout = autoIntensity_Callback(h, eventdata, handles, varargin)
global state gh
genericCallback(h);




% --------------------------------------------------------------------
function varargout = LUTpercent_Callback(h, eventdata, handles, varargin)
global state gh
genericCallback(h);







% --------------------------------------------------------------------
function varargout = highPixelValue_Callback(h, eventdata, handles, varargin)
global gh state
genericCallbackCell(h);
value = get(gh.imageProcessingGUI.fileName, 'Value');
set(state.imageProc.internal.axis{value}, 'CLim', [state.imageProc.lowPixelValue ...
		state.imageProc.highPixelValue]);

if state.imageProc.updateLUTMax 
	set([gca gh.maxProjectionGUI.xyaxis gh.maxProjectionGUI.xzaxis gh.maxProjectionGUI.yzaxis ], 'CLim', [state.imageProc.lowPixelValue ...
			state.imageProc.highPixelValue]);
end


% --------------------------------------------------------------------
function varargout = lowPixelValue_Callback(h, eventdata, handles, varargin)
global gh state
genericCallbackCell(h);
value = get(gh.imageProcessingGUI.fileName, 'Value');
set(state.imageProc.internal.axis{value}, 'CLim', [state.imageProc.lowPixelValue ...
		state.imageProc.highPixelValue]);

if state.imageProc.updateLUTMax 
	set([gca gh.maxProjectionGUI.xyaxis gh.maxProjectionGUI.xzaxis gh.maxProjectionGUI.yzaxis ], 'CLim', [state.imageProc.lowPixelValue ...
			state.imageProc.highPixelValue]);
end



% --------------------------------------------------------------------
function varargout = sumIntensity_Callback(h, eventdata, handles, varargin)

