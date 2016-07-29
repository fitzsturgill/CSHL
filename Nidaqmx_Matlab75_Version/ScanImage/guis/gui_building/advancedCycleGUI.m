function varargout = advancedCycleGUI(varargin)
% ADVANCEDCYCLEGUI M-file for advancedCycleGUI.fig
%      ADVANCEDCYCLEGUI, by itself, creates a new ADVANCEDCYCLEGUI or raises the existing
%      singleton*.
%
%      H = ADVANCEDCYCLEGUI returns the handle to a new ADVANCEDCYCLEGUI or the handle to
%      the existing singleton*.
%
%      ADVANCEDCYCLEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ADVANCEDCYCLEGUI.M with the given input
%      arguments.
%
%      ADVANCEDCYCLEGUI('Property','Value',...) creates a new ADVANCEDCYCLEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before advancedCycleGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to advancedCycleGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help advancedCycleGUI

% Last Modified by GUIDE v2.5 10-Oct-2012 16:17:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @advancedCycleGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @advancedCycleGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before advancedCycleGUI is made visible.
function advancedCycleGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to advancedCycleGUI (see VARARGIN)

% Choose default command line output for advancedCycleGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes advancedCycleGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = advancedCycleGUI_OutputFcn(hObject, eventdata, handles)
	varargout{1} = handles.output;

function deleteCyclePos_Callback(h, eventdata, handles)
	deleteCyclePosition

function insertCyclePos_Callback(h, eventdata, handles)
	insertCyclePosition
	
function generic_Callback(h, eventdata, handles)
%	set(h, 'Enable', 'off');
	genericCallback(h);
	global state
	if state.cycle.displayCyclePosition==state.cycle.currentCyclePosition 
		state.phys.internal.needNewOutputData=1;
	end
	set(h, 'Enable', 'on');
	
function repeatsDone_Callback(h, eventdata, handles)
	genericCallback(h);

function displayCyclePosition_Callback(h, eventdata, handles)
	genericCallback(h);
	global state
	updateCycleDisplay(state.cycle.displayCyclePosition);

function redefineCycle_Callback(h, eventdata, handles)
	genericCallback(h);
%	set(h, 'Enable', 'off');

	try
		global state
		persistent guiOrder
		if isempty(guiOrder) | ~iscell(guiOrder)
			guiOrder={'repeats', 'delay', 'frames', 'blaster', 'imageOn', 'physOn', 'tracker', 'avgFrames', 'linescan', 'da0', 'da1', 'aux4', 'aux5', 'aux6', 'aux7',...
                'mcOn', 'mcOlfOn', 'mcOlfDelay', 'mcOlfDuration', 'mcOlfValve'};
		end
		
		state.internal.cycleChanged=1;
		tag=get(h, 'Tag');
		if ~isempty(findstr(tag, 'Slider'))
			tag=tag(1:end-6);
		end
	
		eval(['state.cycle.' tag 'List(state.cycle.displayCyclePosition)=state.cycle.' tag ';']);
	
		p=strcmp(guiOrder, tag);
		if any(p) 
			global gh
			p=mod(find(p)+1, length(guiOrder));
			if p==0 
				p=length(guiOrder);
			end
			eval(['setfocus(gh.advancedCycleGui.' guiOrder{p} ');']);
		end
	
		if state.cycle.displayCyclePosition==state.cycle.currentCyclePosition 
			applyAdvancedCyclePosition;
		end
	catch
		disp(lasterr);
	end
	set(h, 'Enable', 'on');


% --- Executes on button press in writeProtect.
function writeProtect_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	toggleCycleGUIS;


% --- Executes on button press in randomize.
function randomize_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	setupCycleRandomList;


% --- Executes on button press in useCyclePos.
function useCyclePos_Callback(hObject, eventdata, handles)
	genericCallback(hObject);




% --- Executes on button press in loadCurrentButton.
function loadCurrentButton_Callback(hObject, eventdata, handles)
	putCurrentInCyclePos;





% --- Executes on button press in mcOlfOn.
function mcOlfOn_Callback(hObject, eventdata, handles)
% hObject    handle to mcOlfOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mcOlfOn



function mcOlfDelaySlider_Callback(hObject, eventdata, handles)
% hObject    handle to mcOlfDelaySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mcOlfDelaySlider as text
%        str2double(get(hObject,'String')) returns contents of mcOlfDelaySlider as a double


% --- Executes during object creation, after setting all properties.
function mcOlfDelaySlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfDelaySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider21_Callback(hObject, eventdata, handles)
% hObject    handle to mcOlfDelaySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfDelaySlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% --- Executes during object creation, after setting all properties.
function mcOlfDuration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfDuration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object creation, after setting all properties.
function mcOlfDurationSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfDurationSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% --- Executes during object creation, after setting all properties.
function mcOlfValve_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfValve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function mcOlfValveSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfValveSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function mcOlfDelay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mcOlfDelay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    copyCyclePosition;


% --- Executes on button press in olfShuntEnabled.
function olfShuntEnabled_Callback(hObject, eventdata, handles)
% hObject    handle to olfShuntEnabled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of olfShuntEnabled
    global state
    bool = get(hObject, 'Value');
    genericCallback(hObject);
    if ~bool && length(state.phys.mcAcq.olfDevice.Line) > length(state.phys.mcAcq.olfHwLines)
        delete(state.phys.mcAcq.olfDevice.Line(end));
        state.phys.mcAcq.olfShuntLine=[];
    elseif bool && length(state.phys.mcAcq.olfDevice.Line) == length(state.phys.mcAcq.olfHwLines)
        state.phys.mcAcq.olfShuntLine = addline(state.phys.mcAcq.olfDevice, state.phys.mcAcq.olfShuntHwLine, state.phys.mcAcq.olfPort, 'out');
    end
