function varargout = maxProjectionGUI(varargin)
% MAXPROJECTIONGUI Application M-file for maxProjectionGUI.fig
%    FIG = MAXPROJECTIONGUI launch maxProjectionGUI GUI.
%    MAXPROJECTIONGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Feb-2002 09:25:26

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
function varargout = XYMax_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.
global state gh
set(gh.maxProjectionGUI.XZMax, 'Value', 0);
set(gh.maxProjectionGUI.YZMax, 'Value', 0);
genericCallback(h);
genericCallback(gh.maxProjectionGUI.XZMax);
genericCallback(gh.maxProjectionGUI.YZMax);


% --------------------------------------------------------------------
function varargout = XZMax_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.
global state gh
set(gh.maxProjectionGUI.XYMax, 'Value', 0);
set(gh.maxProjectionGUI.YZMax, 'Value', 0);
genericCallback(h);
genericCallback(gh.maxProjectionGUI.XYMax);
genericCallback(gh.maxProjectionGUI.YZMax);

% --------------------------------------------------------------------
function varargout = YZMax_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.radiobutton3.
global state gh
set(gh.maxProjectionGUI.XYMax, 'Value', 0);
set(gh.maxProjectionGUI.XZMax, 'Value', 0);
genericCallback(h);
genericCallback(gh.maxProjectionGUI.XYMax);
genericCallback(gh.maxProjectionGUI.XZMax);

% --------------------------------------------------------------------
function varargout = generic_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrameSlider.
global state gh
genericCallbackCell(h);

% --------------------------------------------------------------------
function varargout = fileName_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrame.

% --------------------------------------------------------------------
function varargout = maxFileName_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.currentFrame.




% --------------------------------------------------------------------
function varargout = projectall_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.
global state gh 
set(gh.maxProjectionGUI.figure1, 'Pointer', 'watch');
value = get(gh.maxProjectionGUI.fileName, 'Value');

for i = 1:3
    switch i
    case 1
        direction = 'XY';
    case 2
        direction = 'XZ';
    case 3
        direction = 'YZ';
    end
    if state.imageProc.averageNotProject    %average not project
        a{i} = collapse(state.imageProc.cell.currentImage{value}...
            (:,:,state.imageProc.cell.maxStart{value}:state.imageProc.cell.maxEnd{value}), direction,'average');
    else
        a{i} = collapse(state.imageProc.cell.currentImage{value}...
            (:,:,state.imageProc.cell.maxStart{value}:state.imageProc.cell.maxEnd{value}), direction);
    end
    
end

xyimage = image('CData', a{1}, 'CDataMapping', 'scaled', 'parent', gh.maxProjectionGUI.xyaxis);
set(gh.maxProjectionGUI.xyaxis, 'XLim',[1 size(a{1},2)], 'YLim', [1 size(a{1},1)]);
xzimage = image('CData', a{2}, 'CDataMapping', 'scaled', 'parent', gh.maxProjectionGUI.xzaxis);
set(gh.maxProjectionGUI.xzaxis, 'XLim',[1 size(a{2},2)], 'YLim', [1 size(a{2},1)]);
yzimage = image('CData', a{3}, 'CDataMapping', 'scaled', 'parent', gh.maxProjectionGUI.yzaxis);
set(gh.maxProjectionGUI.yzaxis, 'XLim',[1 size(a{3},2)], 'YLim', [1 size(a{3},1)]);

min = min(min(min(state.imageProc.cell.currentImage{value})));
max = .4*double(max(max(max(state.imageProc.cell.currentImage{value}))));
try
    set([gh.maxProjectionGUI.xyaxis gh.maxProjectionGUI.xzaxis gh.maxProjectionGUI.yzaxis], 'Clim', [min max]);
end

colormap(makeColorMap('gray',8));
set(gh.maxProjectionGUI.figure1, 'Pointer', 'arrow');
state.imageProc.internal.maxValue = value;
state.imageProc.internal.maxProj = a;

% --------------------------------------------------------------------
function varargout = threeDCube_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.pushbutton2.
global state gh
set(gh.maxProjectionGUI.figure1, 'Pointer', 'watch');
value = get(gh.maxProjectionGUI.fileName, 'Value');

if value ~= state.imageProc.internal.maxValue
    
    for i = 1:3
        switch i
        case 1
            direction = 'XY';
        case 2
            direction = 'XZ';
        case 3
            direction = 'YZ';
        end
        if state.imageProc.averageNotProject    %average not project
            a{i} = collapse(state.imageProc.cell.currentImage{value}...
                (:,:,state.imageProc.cell.maxStart{value}:state.imageProc.cell.maxEnd{value}), direction,'average');
        else
            a{i} = collapse(state.imageProc.cell.currentImage{value}...
                (:,:,state.imageProc.cell.maxStart{value}:state.imageProc.cell.maxEnd{value}), direction);
        end
        
    end
    xyimage = image('CData', a{1}, 'CDataMapping', 'scaled', 'parent', gh.maxProjectionGUI.xyaxis);
    set(gh.maxProjectionGUI.xyaxis, 'XLim',[1 size(a{1},2)], 'YLim', [1 size(a{1},1)]);
    xzimage = image('CData', a{2}, 'CDataMapping', 'scaled', 'parent', gh.maxProjectionGUI.xzaxis);
    set(gh.maxProjectionGUI.xzaxis, 'XLim',[1 size(a{2},2)], 'YLim', [1 size(a{2},1)]);
    yzimage = image('CData', a{3}, 'CDataMapping', 'scaled', 'parent', gh.maxProjectionGUI.yzaxis);
    set(gh.maxProjectionGUI.yzaxis, 'XLim',[1 size(a{3},2)], 'YLim', [1 size(a{3},1)]);
    
    min = min(min(min(state.imageProc.cell.currentImage{value})));
    max = .4*double(max(max(max(state.imageProc.cell.currentImage{value}))));
    try
        set([gh.maxProjectionGUI.xyaxis gh.maxProjectionGUI.xzaxis gh.maxProjectionGUI.yzaxis], 'Clim', [min max]);
    end
    colormap(makeColorMap('gray',8));	
else
    a = state.imageProc.internal.maxProj;
end

frames = size(a{2},1);
rows = size(a{1},1);
columns = size(a{1},2);
if frames > 100 % binn before displaying
    a{1}=genericBin(a{1},rows/128,columns/128,1);
    a{2}=genericBin(a{2},rows/128,1,1);
    a{3}=genericBin(a{3},1,columns/128,1);
    rows = size(a{1},1);
    columns = size(a{1},2);
end

b = cat(3,a{1},zeros(rows,columns,(frames-2)),a{1});
b(1,:,:) = flipdim(a{2}',2);
b(1,:,:) = flipdim(b(1,:,:),1);
b(end,:,:) = b(1,:,:);
b(:,1,:) = a{3};
b(:,1,:) = flipdim(b(:,1,:),3);
b(:,end,:) = b(:,1,:);
b = double(b);
state.imageProc.cubeHandle= cubeit(b,rows,columns,frames);
set(gh.maxProjectionGUI.figure1, 'Pointer', 'arrow');


% --------------------------------------------------------------------
function varargout = maxRatioY_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);
try
    set(state.imageProc.cubeaxis, 'DataAspectRatio', [state.imageProc.maxRatioX1 state.imageProc.maxRatioY state.imageProc.maxRatioZ]);
end
try 
    set(state.imageProc.threedaxis, 'DataAspectRatio', [state.imageProc.maxRatioX1 state.imageProc.maxRatioY state.imageProc.maxRatioZ]);
end
% --------------------------------------------------------------------
function varargout = maxRatioZ_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);
try
    set(state.imageProc.cubeaxis, 'DataAspectRatio', [state.imageProc.maxRatioX1 state.imageProc.maxRatioY state.imageProc.maxRatioZ]);
end
try 
    set(state.imageProc.threedaxis, 'DataAspectRatio', [state.imageProc.maxRatioX1 state.imageProc.maxRatioY state.imageProc.maxRatioZ]);
end

% --------------------------------------------------------------------
function varargout = updateLUTMax_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);


% --------------------------------------------------------------------
function varargout = maxRatioX1_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);
try 
    set(state.imageProc.cubeaxis, 'DataAspectRatio', [state.imageProc.maxRatioX1 state.imageProc.maxRatioY state.imageProc.maxRatioZ]);
end
try 
    set(state.imageProc.threedaxis, 'DataAspectRatio', [state.imageProc.maxRatioX1 state.imageProc.maxRatioY state.imageProc.maxRatioZ]);
end



% --------------------------------------------------------------------
function varargout = autoShrink_Callback(h, eventdata, handles, varargin)
global gh state
genericCallback(h);


% --------------------------------------------------------------------
function varargout = threeDRecon_Callback(h, eventdata, handles, varargin)
global gh state
do3DReconstruction;



% --------------------------------------------------------------------
function varargout = pushbutton2_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = averageNotProject_Callback(h, eventdata, handles, varargin)
genericCallback(h);

