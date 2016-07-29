function varargout = spineGUI(varargin)
% SPINEGUI Application M-file for spineGUI.fig
%    FIG = SPINEGUI launch spineGUI GUI.
%    SPINEGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 09-May-2002 11:53:01

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
function varargout = genericMax_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallback(h);
state.imageProc.spine.startMax = round(state.imageProc.spine.startMax);
updateGUIByGlobal('state.imageProc.spine.startMax');
state.imageProc.spine.stopMax = round(state.imageProc.spine.stopMax);
updateGUIByGlobal('state.imageProc.spine.stopMax');

function varargout = genericAuto_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallback(h);


% --------------------------------------------------------------------
function varargout = generic_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallback(h);
state.imageProc.spine.lowPixelValue=round(double(state.imageProc.spine.lowPixelValue));
state.imageProc.spine.highPixelValue=round(double(state.imageProc.spine.highPixelValue));
updateGuiByGlobal('state.imageProc.spine.highPixelValue');
updateGuiByGlobal('state.imageProc.spine.lowPixelValue');

if state.imageProc.spine.topImage
	set([gh.spineGUI.mainAxes gh.spineGUI.initialaxis], 'Clim', [state.imageProc.spine.lowPixelValue ...
			state.imageProc.spine.highPixelValue]);
else
	set(gh.spineGUI.initialaxis2, 'Clim', [state.imageProc.spine.lowPixelValue ...
			state.imageProc.spine.highPixelValue]);
end


% --------------------------------------------------------------------
function varargout = genericThreshold_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallback(h);

% --------------------------------------------------------------------
function varargout = genericRatio_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallback(h);




% --------------------------------------------------------------------
function varargout = currentSpineSlider_Callback(h, eventdata, handles, varargin)
global gh state
% Stub for Callback of the uicontrol handles.currentSpineSlider.

state.imageProc.spine.maxFlag=0;
if state.imageProc.spine.maxFlag == 1
	return
end
genericCallback(h);
state.imageProc.spine.currentSpineFrame = round(state.imageProc.spine.currentSpineFrame);
updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
set(state.imageProc.spine.imagehandle, 'CData', state.imageProc.spine.initialImage(:,:,state.imageProc.spine.currentSpineFrame));
state.imageProc.spine.ZCurrentTop = state.imageProc.spine.ZStartTop + (state.imageProc.spine.currentSpineFrame-1)*state.imageProc.spine.zStepSizeTop;
updateGUIByGlobal('state.imageProc.spine.ZCurrentTop');

if state.imageProc.spine.autoTransferTop
	if state.imageProc.spine.auto
		showBWImage;
		return
	end
	
	if state.imageProc.spine.eraseSpines
		try
			removeLines;
		end
	end
	state.imageProc.spine.topImage=1;
	updateGuiByGlobal('state.imageProc.spine.topImage');

	state.imageProc.spine.bottomImage=1;
	updateGuiByGlobal('state.imageProc.spine.bottomImage');
	transfer_Callback;
end

% --------------------------------------------------------------------

function varargout = currentSpineSliderBot_Callback(h, eventdata, handles, varargin)
global gh state
state.imageProc.spine.maxFlag2=0;
if state.imageProc.spine.maxFlag2 == 1
	return
end
genericCallback(h);
state.imageProc.spine.currentSpineFrame2 = round(state.imageProc.spine.currentSpineFrame2);
updateGUIByGlobal('state.imageProc.spine.currentSpineFrame2');
set(state.imageProc.spine.imagehandle2, 'CData', state.imageProc.spine.initialImage2(:,:,state.imageProc.spine.currentSpineFrame2));
state.imageProc.spine.ZCurrentBot = state.imageProc.spine.ZStartBot + (state.imageProc.spine.currentSpineFrame2-1)*state.imageProc.spine.zStepSizeBot;
updateGUIByGlobal('state.imageProc.spine.ZCurrentBot');

if state.imageProc.spine.autoTransferBot
	if state.imageProc.spine.auto
		showBWImage;
		return
	end
	
	if state.imageProc.spine.eraseSpines
		try
			removeLines;
		end
	end
	
	state.imageProc.spine.baseName=[state.imageProc.spine.loadedFileNameBot '0'];
	updateGUIByGlobal('state.imageProc.spine.baseName');
	
	if state.imageProc.spine.maxFlag2 == 0 % max not dsipalyed
		[x,x1,y,y1] = getCurrentAxisLimits(gh.spineGUI.initialaxis2);
		set(gh.spineGUI.mainAxes, 'YLim', [y y1], 'XLim', [x x1], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
				state.imageProc.spine.highPixelValue]);
		
		state.imageProc.spine.mainImage = image('CData', state.imageProc.spine.initialImage2(y:y1,x:x1,state.imageProc.spine.currentSpineFrame2),...
			'CDataMapping', 'scaled','Parent', gh.spineGUI.mainAxes,'YData', [y y1], 'XData', [x x1]);
		state.imageProc.spine.from = state.imageProc.spine.ZCurrentBot;
		state.imageProc.spine.to = state.imageProc.spine.ZCurrentBot;
	else
		[x,x1,y,y1] = getCurrentAxisLimits(gh.spineGUI.initialaxis2);
		set(gh.spineGUI.mainAxes, 'YLim', [y y1], 'XLim', [x x1], 'Ydir', 'reverse', 'Clim', [state.imageProc.spine.lowPixelValue ...
				state.imageProc.spine.highPixelValue]);
		
		state.imageProc.spine.mainImage = image('CData', state.imageProc.spine.maxProjection2(y:y1,x:x1),...
			'CDataMapping', 'scaled','Parent', gh.spineGUI.mainAxes,'YData', [y y1], 'XData', [x x1]);
		state.imageProc.spine.from = state.imageProc.spine.from2;
		state.imageProc.spine.to = state.imageProc.spine.to2;
	end
	reshuffleAxisHandles(gh.spineGUI.mainAxes);
end





% --------------------------------------------------------------------
function varargout = transfer_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.transferCurrentImage.
	transferImages

% --------------------------------------------------------------------
function varargout = preview_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.transferCurrentImage.
global gh state
% Call load image ui
loadPreviewSpine;

% --------------------------------------------------------------------
function varargout = project_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.transferCurrentImage.
global gh state
set(gh.spineGUI.figure1, 'Pointer', 'watch');


if get(gh.spineGUI.XY, 'Value')
	direction = get(gh.spineGUI.XY, 'String');
elseif get(gh.spineGUI.XZ, 'Value')
	direction = get(gh.spineGUI.XZ, 'String');
elseif get(gh.spineGUI.YZ, 'Value')
	direction = get(gh.spineGUI.YZ, 'String');
else
	display('No Direction Specified: Check state.imageProc.internal.');
	return	
end

if state.imageProc.spine.topImage
% 	if state.imageProc.spine.maxFlag == 1
% 		reloadSpineImage;
% 	end
	state.imageProc.spine.reloadFrameTop = state.imageProc.spine.currentSpineFrame;
	state.imageProc.spine.from1 = state.imageProc.spine.ZStartTop + (state.imageProc.spine.startMax-1)*state.imageProc.spine.zStepSizeTop;
	state.imageProc.spine.to1 = state.imageProc.spine.ZStartTop + (state.imageProc.spine.stopMax-1)*state.imageProc.spine.zStepSizeTop;
	state.imageProc.spine.maxProjection = collapse(state.imageProc.spine.initialImage(:,:,state.imageProc.spine.startMax:state.imageProc.spine.stopMax), direction);
	state.imageProc.spine.maxFlag = 1;
	columns = size(state.imageProc.spine.maxProjection,2);
	rows = size(state.imageProc.spine.maxProjection,1);
% 	state.imageProc.spine.stopMax = 1.0001;
% 	state.imageProc.spine.totalSpineFrames = state.imageProc.spine.stopMax;
% 	state.imageProc.spine.currentSpineFrame = 1;
% 	state.imageProc.spine.startMax = 1;
% 	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
% 	updateGUIByGlobal('state.imageProc.spine.stopMax');
% 	updateGUIByGlobal('state.imageProc.spine.startMax');
% 	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames');
	
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', (state.imageProc.spine.totalSpineFrames+.001), 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	
% 	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.maxProjection)));
% 	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
% 	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.maxProjection))));
% 	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	if strcmp(direction,'XY') 
		set(gh.spineGUI.initialaxis,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
			[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	else
		set(gh.spineGUI.initialaxis, 'Ydir', 'reverse', 'CLim', ...
			[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	end
	set(state.imageProc.spine.imagehandle,'CData', state.imageProc.spine.maxProjection(:,:,1), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis, 'ButtonDownFcn', 'spineImageOverFcn');
	
	if state.imageProc.spine.autoTransferTop
		if state.imageProc.spine.auto
			showBWImage;
			return
		end
		
		if state.imageProc.spine.eraseSpines
			try
				removeLines;
			end
		end
		state.imageProc.spine.topImage=1;
		updateGuiByGlobal('state.imageProc.spine.topImage');
	
		state.imageProc.spine.bottomImage=1;
		updateGuiByGlobal('state.imageProc.spine.bottomImage');
		transfer_Callback;
	end
	
elseif state.imageProc.spine.bottomImage
	if state.imageProc.spine.maxFlag2 == 1
		reloadSpineImage;
	end
	state.imageProc.spine.reloadFrameBot = state.imageProc.spine.currentSpineFrame2;
	state.imageProc.spine.from2 = state.imageProc.spine.ZStartBot + (state.imageProc.spine.startMax-1)*state.imageProc.spine.zStepSizeBot;
	state.imageProc.spine.to2 = state.imageProc.spine.ZStartBot + (state.imageProc.spine.stopMax-1)*state.imageProc.spine.zStepSizeBot;
	state.imageProc.spine.maxProjection2 = collapse(state.imageProc.spine.initialImage2(:,:,state.imageProc.spine.startMax:state.imageProc.spine.stopMax), direction);
	state.imageProc.spine.maxFlag2 = 1;
	columns = size(state.imageProc.spine.maxProjection2,2);
	rows = size(state.imageProc.spine.maxProjection2,1);
	state.imageProc.spine.stopMax = 1.0001;
	state.imageProc.spine.totalSpineFrames2 = state.imageProc.spine.stopMax;
	state.imageProc.spine.currentSpineFrame2 = 1;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame2');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames2');
	
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider2], 'Min', 1, 'Max', (state.imageProc.spine.totalSpineFrames2+.001), 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	
	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.maxProjection2)));
	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.maxProjection2))));
	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	if strcmp(direction,'XY') 
		set(gh.spineGUI.initialaxis2,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
			[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	else
		set(gh.spineGUI.initialaxis2, 'Ydir', 'reverse', 'CLim', ...
			[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	end
	
	set(state.imageProc.spine.imagehandle2, 'CData', state.imageProc.spine.maxProjection2(:,:,1), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis2, 'ButtonDownFcn', 'spineImageOverFcn');
end

set(gh.spineGUI.figure1, 'Pointer', 'arrow');

function varargout = reload_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.XY.
global gh state
reloadSpineImage;


% --------------------------------------------------------------------
function varargout = XY_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.XY.
global gh state
set(gh.spineGUI.XZ, 'Value', 0);
set(gh.spineGUI.YZ, 'Value', 0);
genericCallback(h);
genericCallback(gh.spineGUI.XZ);
genericCallback(gh.spineGUI.YZ);

% --------------------------------------------------------------------
function varargout = XZ_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.XZ.
global gh state
set(gh.spineGUI.YZ, 'Value', 0);
set(gh.spineGUI.XY, 'Value', 0);
genericCallback(h);
genericCallback(gh.spineGUI.YZ);
genericCallback(gh.spineGUI.XY);

% --------------------------------------------------------------------
function varargout = YZ_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.YZ.
global gh state
set(gh.spineGUI.XZ, 'Value', 0);
set(gh.spineGUI.XY, 'Value', 0);
genericCallback(h);
genericCallback(gh.spineGUI.XZ);
genericCallback(gh.spineGUI.XY);

% --------------------------------------------------------------------
function varargout = top_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.YZ.
global gh state
set(gh.spineGUI.bottomImage, 'Value', 0);
genericCallback(h);
genericCallback(gh.spineGUI.bottomImage);

if state.imageProc.spine.topImage
	clim = get(gh.spineGUI.initialaxis, 'CLim');
	state.imageProc.spine.lowPixelValue = clim(1);
	state.imageProc.spine.highPixelValue = clim(2);
	updateGUIByGlobAL('state.imageProc.spine.highPixelValue');
	updateGUIByGlobAL('state.imageProc.spine.lowPixelValue');
	columns = size(state.imageProc.spine.initialImage,2);
	rows = size(state.imageProc.spine.initialImage,1);
	if state.imageProc.spine.maxFlag == 0
		state.imageProc.spine.totalSpineFrames = size(state.imageProc.spine.initialImage,3);
	else
		state.imageProc.spine.totalSpineFrames =1.001;
	end
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames');
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider ], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	set(gh.spineGUI.initialaxis, 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	
	
elseif state.imageProc.spine.bottomImage
	clim = get(gh.spineGUI.initialaxis2, 'CLim');
	state.imageProc.spine.lowPixelValue = clim(1);
	state.imageProc.spine.highPixelValue = clim(2);
	updateGUIByGlobAL('state.imageProc.spine.highPixelValue');
	updateGUIByGlobAL('state.imageProc.spine.lowPixelValue');
	columns = size(state.imageProc.spine.initialImage2,2);
	rows = size(state.imageProc.spine.initialImage2,1);
	if state.imageProc.spine.maxFlag2 == 0
		state.imageProc.spine.totalSpineFrames2 = size(state.imageProc.spine.initialImage2,3);
	else
		state.imageProc.spine.totalSpineFrames2 =1.001;
	end
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames2;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames2');
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames2, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	set(gh.spineGUI.initialaxis2, 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);	
else
	display('No Image Loaded');
	return
end
% --------------------------------------------------------------------
function varargout = bottom_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.YZ.
global gh state
set(gh.spineGUI.topImage, 'Value', 0);
genericCallback(h);
genericCallback(gh.spineGUI.topImage);

if state.imageProc.spine.topImage
	clim = get(gh.spineGUI.initialaxis, 'CLim');
	state.imageProc.spine.lowPixelValue = clim(1);
	state.imageProc.spine.highPixelValue = clim(2);
	updateGUIByGlobAL('state.imageProc.spine.highPixelValue');
	updateGUIByGlobAL('state.imageProc.spine.lowPixelValue');
	columns = size(state.imageProc.spine.initialImage,2);
	rows = size(state.imageProc.spine.initialImage,1);
	if state.imageProc.spine.maxFlag == 0
		state.imageProc.spine.totalSpineFrames = size(state.imageProc.spine.initialImage,3);
	else
		state.imageProc.spine.totalSpineFrames =1.001;
	end
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames');
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider ], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	set(gh.spineGUI.initialaxis, 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	
	
elseif state.imageProc.spine.bottomImage
	clim = get(gh.spineGUI.initialaxis2, 'CLim');
	state.imageProc.spine.lowPixelValue = clim(1);
	state.imageProc.spine.highPixelValue = clim(2);
	updateGUIByGlobAL('state.imageProc.spine.highPixelValue');
	updateGUIByGlobAL('state.imageProc.spine.lowPixelValue');
	columns = size(state.imageProc.spine.initialImage2,2);
	rows = size(state.imageProc.spine.initialImage2,1);
	if state.imageProc.spine.maxFlag2 == 0
		state.imageProc.spine.totalSpineFrames2 = size(state.imageProc.spine.initialImage2,3);
	else
		state.imageProc.spine.totalSpineFrames2 =1.001;
	end
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames2;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames2');
	set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames2, 'SliderStep', ...
		[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	set(gh.spineGUI.initialaxis2, 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);	
else
	display('No Image Loaded');
	return
end
% --------------------------------------------------------------------
function varargout = filter_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.auto.

global gh state
set(gh.spineGUI.figure1, 'Pointer', 'watch');

if state.imageProc.spine.topImage
	
	sizeArray = size(state.imageProc.spine.initialImage,3);
	h = waitbar(0,'Filtering Image...', 'Name', 'Filter Status', 'Pointer', 'watch');
	
	for i = 1:sizeArray
		filteredImage(:,:,i) = medfilt2(state.imageProc.spine.initialImage(:,:,i),[3 3]);
		waitbar(i/sizeArray,h, ['Filtering Frame Number ' num2str(i)]);
	end
	
	state.imageProc.spine.maxFlag = 0;
	state.imageProc.spine.initialImage = filteredImage;
	columns = size(state.imageProc.spine.initialImage,2);
	rows = size(state.imageProc.spine.initialImage,1);
	state.imageProc.spine.totalSpineFrames = size(state.imageProc.spine.initialImage,3);
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames;
	state.imageProc.spine.currentSpineFrame = 1;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames');
	if state.imageProc.spine.totalSpineFrames == 1
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', 1.001, 'SliderStep', ...
			[1 1]);
	else
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames, 'SliderStep', ...
			[1/state.imageProc.spine.totalSpineFrames 1/state.imageProc.spine.totalSpineFrames]);
	end
	
	
% 	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.initialImage)));
% 	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
% 	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.initialImage))));
% 	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	set(gh.spineGUI.initialaxis,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	state.imageProc.spine.XLim = [1 columns];
	state.imageProc.spine.YLim = [1 rows];
	state.imageProc.spine.imagehandle = image('CData', state.imageProc.spine.initialImage(:,:,1), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis, 'ButtonDownFcn', 'spineImageOverFcn');
	close(h);
	
elseif state.imageProc.spine.bottomImage
	
	sizeArray = size(state.imageProc.spine.initialImage2,3);
	h = waitbar(0,'Filtering Image...', 'Name', 'Filter Status', 'Pointer', 'watch');
	
	for i = 1:sizeArray
		filteredImage(:,:,i) = medfilt2(state.imageProc.spine.initialImage2(:,:,i),[3 3]);
		waitbar(i/sizeArray,h, ['Filtering Frame Number ' num2str(i)]);
	end
	
	state.imageProc.spine.maxFlag = 0;
	state.imageProc.spine.initialImage2 = filteredImage;
	columns = size(state.imageProc.spine.initialImage2,2);
	rows = size(state.imageProc.spine.initialImage2,1);
	state.imageProc.spine.totalSpineFrames2 = size(state.imageProc.spine.initialImage2,3);
	state.imageProc.spine.stopMax = state.imageProc.spine.totalSpineFrames2;
	state.imageProc.spine.currentSpineFrame2 = 1;
	state.imageProc.spine.startMax = 1;
	updateGUIByGlobal('state.imageProc.spine.currentSpineFrame2');
	updateGUIByGlobal('state.imageProc.spine.startMax');
	updateGUIByGlobal('state.imageProc.spine.stopMax');
	updateGUIByGlobal('state.imageProc.spine.totalSpineFrames2');
	
	if state.imageProc.spine.totalSpineFrames2 == 1
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', 1.001, 'SliderStep', ...
			[1 1]);
	else
		set([gh.spineGUI.startMaxSlider gh.spineGUI.stopMaxSlider gh.spineGUI.currentSpineSlider], 'Min', 1, 'Max', state.imageProc.spine.totalSpineFrames, 'SliderStep', ...
			[1/state.imageProc.spine.totalSpineFrames2 1/state.imageProc.spine.totalSpineFrames2]);
	end
	state.imageProc.spine.lowPixelValue = min(min(min(state.imageProc.spine.initialImage2)));
	updateGUIByGlobal('state.imageProc.spine.lowPixelValue');
	state.imageProc.spine.highPixelValue = .1*double(max(max(max(state.imageProc.spine.initialImage2))));
	updateGUIByGlobal('state.imageProc.spine.highPixelValue');
	
	set(gh.spineGUI.figure1,'Colormap', makeColormap('gray',8));
	set(gh.spineGUI.initialaxis2,'YLim', [1 rows], 'XLim', [1 columns], 'Ydir', 'reverse', 'CLim', ...
		[state.imageProc.spine.lowPixelValue state.imageProc.spine.highPixelValue]);
	state.imageProc.spine.XLim = [1 columns];
	state.imageProc.spine.YLim = [1 rows];
	state.imageProc.spine.imagehandle2 = image('CData', state.imageProc.spine.initialImage2(:,:,1), ...
		'CDataMapping', 'scaled', 'Parent', gh.spineGUI.initialaxis2, 'ButtonDownFcn', 'spineImageOverFcnBot');
	close(h);
else
	return
end
set(gh.spineGUI.figure1, 'Pointer', 'arrow');

% --------------------------------------------------------------------
function varargout = auto_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.auto.
global gh state
genericCallback(h);
state.imageProc.spine.manual = 0;
updateGUIByGlobal('state.imageProc.spine.manual');
switchManualAuto;

% --------------------------------------------------------------------
function varargout = manual_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.manual.
global gh state
genericCallback(h);
state.imageProc.spine.auto = 0;
updateGUIByGlobal('state.imageProc.spine.auto');
switchManualAuto;


% --------------------------------------------------------------------
function varargout = zoomIn_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.zoomIn.
global gh state
zoomIn(gh.spineGUI.initialaxis);

% --------------------------------------------------------------------
function varargout = zoomOut_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.zoomOut.
global gh state
zoomOut(gh.spineGUI.initialaxis);
if state.imageProc.spine.autoTransferTop
	transferImages;
end

% --------------------------------------------------------------------
function varargout = zoomIn2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.zoomIn.
global gh state

zoomIn(gh.spineGUI.initialaxis2);
if state.imageProc.spine.autoTransferTop
	transferImages;
end

% --------------------------------------------------------------------
function varargout = zoomOut2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.zoomOut.
global gh state
set(gh.spineGUI.initialaxis2, 'XLim', state.imageProc.spine.XLim2, 'YLim', state.imageProc.spine.YLim2);
% --------------------------------------------------------------------
function varargout = zoomInPreview_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.zoomIn.
global gh state
zoomIn(gh.spineGUI.previewaxis);

% --------------------------------------------------------------------
function varargout = zoomOutPreview_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.zoomOut.
global gh state
set(gh.spineGUI.previewaxis, 'XLim', state.imageProc.spine.XLim3, 'YLim', state.imageProc.spine.YLim3);




% --------------------------------------------------------------------
function varargout = copy1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.copy1.
global gh
copyToClip(gh.spineGUI.initFigure);


% --------------------------------------------------------------------
function varargout = copy2_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.copy2.
global gh
copyToClip(gh.spineGUI.initFigure2);


% --------------------------------------------------------------------
function varargout = copyMain_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.copyMain.
global gh
copyToClip(gh.spineGUI.mainFigure);


% --------------------------------------------------------------------
function varargout = copyPreview_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.copyPreview.
global gh
copyToClip(gh.spineGUI.previewFigure);


% --------------------------------------------------------------------
function varargout = autoTransferTop_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.autoTransferTop.
genericCallback(h);


% --------------------------------------------------------------------
function varargout = autoTransferBot_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.autoTransferBot.
genericCallback(h);






% --------------------------------------------------------------------
function varargout = reloadButton_Callback(h, eventdata, handles, varargin)
reloadSpineImage;	






% --------------------------------------------------------------------
function varargout = eraseSpines_Callback(h, eventdata, handles, varargin)
genericCallback(h);




% --------------------------------------------------------------------
function varargout = spineLine_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global state gh
try
	set(state.imageProc.spine.linehandles,'linewidth', state.imageProc.spine.spineLine);
end
reshuffleAxisHandles(gh.spineGUI.mainAxes);


% --------------------------------------------------------------------
function varargout = spineText_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global state gh
try
	set(state.imageProc.spine.texthandles,'fontsize', state.imageProc.spine.spineText);
	set(state.imageProc.spine.linehandles,'linewidth', state.imageProc.spine.spineLine);
end
reshuffleAxisHandles(gh.spineGUI.mainAxes);


% --------------------------------------------------------------------
function varargout = spineColorMap_Callback(h, eventdata, handles, varargin)
global state gh

val=get(gh.spineGUI.spineColorMap, 'Value');
str=get(gh.spineGUI.spineColorMap, 'String');
str = str{val};
try
	numberOfSpines=length(state.imageProc.spine.texthandles);
	switch str
	case 'Jet'
		map = jet(numberOfSpines);
	case 'Green'
		map = repmat([0 1 0],numberOfSpines,1);
	case 'Blue'
		map = repmat([0 0 1],numberOfSpines,1);
	case 'Red '
		map = repmat([1 0 0],numberOfSpines,1);
	case 'Gray'
		map = gray(numberOfSpines);
	otherwise
		map = jet(numberOfSpines);
	end
	
	for i = 1:length(state.imageProc.spine.texthandles)
		set(state.imageProc.spine.texthandles(i),'color', map(i,:));
		set(state.imageProc.spine.linehandles(i),'color', map(i,:));
	end
end


% --------------------------------------------------------------------
function varargout = calcSpineVols_Callback(h, eventdata, handles, varargin)
genericCallback(h);
