function varargout = mcAcquisition(varargin)
% MCACQUISITION M-file for mcAcquisition.fig
%      MCACQUISITION, by itself, creates a new MCACQUISITION or raises the existing
%      singleton*.
%
%      H = MCACQUISITION returns the handle to a new MCACQUISITION or the handle to
%      the existing singleton*.
%
%      MCACQUISITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MCACQUISITION.M with the given input arguments.
%
%      MCACQUISITION('Property','Value',...) creates a new MCACQUISITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mcAcquisition_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mcAcquisition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mcAcquisition

% Last Modified by GUIDE v2.5 10-Oct-2012 16:15:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mcAcquisition_OpeningFcn, ...
                   'gui_OutputFcn',  @mcAcquisition_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mcAcquisition is made visible.
function mcAcquisition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mcAcquisition (see VARARGIN)

% Choose default command line output for mcAcquisition
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mcAcquisition wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mcAcquisition_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes during object creation, after setting all properties.
function mcInputRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcInputRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function mcInputRate_Callback(hObject, eventdata, handles)
% hObject    handle to mcInputRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mcInputRate as text
%        str2double(get(hObject,'String')) returns contents of mcInputRate as a double
    genericCallback(hObject);



function mcNChannels_Callback(hObject, eventdata, handles)
% hObject    handle to mcNChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mcNChannels as text
%        str2double(get(hObject,'String')) returns contents of mcNChannels as a double
    genericCallback(hObject);

% --- Executes during object creation, after setting all properties.
function mcNChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcNChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentChannel_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentChannel as text
%        str2double(get(hObject,'String')) returns contents of currentChannel as a double
    global state
    val = str2num(get(hObject, 'String'));
    if val > state.phys.mcAcq.totalChannels
        val = state.phys.mcAcq.totalChannels;
    end    
    mcAcqFlipCurrentChannel(val);




% --- Executes during object creation, after setting all properties.
function currentChannel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentChannel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function currentChannelSlider_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    global state
%     if isempty(state.phys.mcAcq.   % what should this condition be?
%             return
%     end
    val = get(hObject, 'Value');
    mcAcqFlipCurrentChannel(val);

% --- Executes during object creation, after setting all properties.
function currentChannelSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentChannelSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function currentChannelName_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentChannelName as text
%        str2double(get(hObject,'String')) returns contents of currentChannelName as a double
    genericCallback(hObject);
    mcAcqUpdateChannelStruct;

% --- Executes during object creation, after setting all properties.
function currentChannelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentChannelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in currentChannelShow.
function currentChannelShow_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of currentChannelShow
    genericCallback(hObject);
    mcAcqUpdateChannelStruct;

% --- Executes on button press in currentChannelShowFilter.
function currentChannelShowFilter_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelShowFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of currentChannelShowFilter
    genericCallback(hObject);
    mcAcqUpdateChannelStruct;


function currentChannelLowPass_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelLowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentChannelLowPass as text
%        str2double(get(hObject,'String')) returns contents of currentChannelLowPass as a double
    genericCallback(hObject);
    mcAcqUpdateChannelStruct;

% --- Executes during object creation, after setting all properties.
function currentChannelLowPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentChannelLowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function currentChannelHighPass_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of currentChannelHighPass as text
%        str2double(get(hObject,'String')) returns contents of currentChannelHighPass as a double
    genericCallback(hObject);
    mcAcqUpdateChannelStruct;

% --- Executes during object creation, after setting all properties.
function currentChannelHighPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to currentChannelHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mcAcqMakeChannelPlots;


function totalChannels_Callback(hObject, eventdata, handles)
% hObject    handle to totalChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of totalChannels as text
%        str2double(get(hObject,'String')) returns contents of totalChannels as a double
    genericCallback(hObject);


% --- Executes during object creation, after setting all properties.
function totalChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to totalChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Probes_Callback(hObject, eventdata, handles)
% hObject    handle to Probes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    mcAcqCurrentChannelSetAll;



function globalLowPass_Callback(hObject, eventdata, handles)
% hObject    handle to globalLowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of globalLowPass as text
%        str2double(get(hObject,'String')) returns contents of globalLowPass as a double
    genericCallback(hObject);


% --- Executes during object creation, after setting all properties.
function globalLowPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to globalLowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function globalHighPass_Callback(hObject, eventdata, handles)
% hObject    handle to globalHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of globalHighPass as text
%        str2double(get(hObject,'String')) returns contents of globalHighPass as a double
    genericCallback(hObject);


% --- Executes during object creation, after setting all properties.
function globalHighPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to globalHighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in currentChannelMUAInclude.
function currentChannelMUAInclude_Callback(hObject, eventdata, handles)
% hObject    handle to currentChannelMUAInclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of currentChannelMUAInclude
    global state
    genericCallback(hObject);
    mcAcqUpdateChannelStruct;
    if ~state.phys.mcAcq.MUA.resetMUA
        state.phys.mcAcq.MUA.resetMUA=1;
        updateGUIByGlobal('state.phys.mcAcq.MUA.resetMUA');
        mcAcqResetMUA;
    end





% --- Executes on button press in resetMUA.
function resetMUA_Callback(hObject, eventdata, handles)
% hObject    handle to resetMUA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of resetMUA
    genericCallback(hObject);
    mcAcqResetMUA;


% --- Executes on button press in olfShuntEnabled.
function olfShuntEnabled_Callback(hObject, eventdata, handles)
% hObject    handle to olfShuntEnabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of olfShuntEnabled
    genericCallback(hObject);
