function varargout = blaster(varargin)
% BLASTER M-file for blaster.fig
%      BLASTER, by itself, creates a new BLASTER or raises the existing
%      singleton*.
%
%      H = BLASTER returns the handle to a new BLASTER or the handle to
%      the existing singleton*.
%
%      BLASTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BLASTER.M with the given input arguments.
%
%      BLASTER('Property','Value',...) creates a new BLASTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before blaster_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to blaster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help blaster

% Last Modified by GUIDE v2.5 13-Dec-2005 14:36:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @blaster_OpeningFcn, ...
                   'gui_OutputFcn',  @blaster_OutputFcn, ...
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


% --- Executes just before blaster is made visible.
function blaster_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = blaster_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function varargout = blankImaging_Callback(hObject, eventdata, handles, varargin)
	set(hObject, 'Enable', 'off');	
	try
		genericCallback(hObject);
		global state
		state.internal.needNewPcellPowerOutput=1;
		handleChange;
	catch
		disp(['blaster: ' lasterr]);
	end
	set(hObject, 'Enable', 'on');	

function active_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	set(hObject, 'Enable', 'off');	
	try
		global state
		if state.blaster.blankImaging
			state.internal.needNewPcellPowerOutput=1;
		end

		state.internal.needNewPcellRepeatedOutput=1;
		state.internal.needNewRepeatedMirrorOutput=1;
		applyChangesToOutput;
	catch
		disp(['blaster: ' lasterr]);
	end
	try
		updateReferenceImage;
	catch
		disp(['blaster (updateReferenceImage): ' lasterr]);
	end
	set(hObject, 'Enable', 'on');	

function varargout = genericConfig_Callback(hObject, eventdata, handles, varargin)
	set(hObject, 'Enable', 'off');	
	try
		genericCallback(hObject);
		global state

		% is the blaster on and the position in use in the current config?
		needUpdate=state.blaster.active & (state.blaster.currentConfig == state.blaster.displayConfig);
			%& any(state.blaster.displayPos==state.blaster.allConfigs{state.blaster.currentConfig, 2}(:,1));	

		switch get(hObject, 'Tag')
			case {'linePos', 'linePosSlider'}		% position 1 in the config
				state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 1)=state.blaster.linePos;
				if needUpdate
					state.internal.needNewRepeatedMirrorOutput=1;
				end
			case {'linePat', 'linePatSlider'}		% position 2 in the config
				state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 2)=state.blaster.linePat;
				if needUpdate
					state.internal.needNewPcellRepeatedOutput=1;
					state.internal.needNewRepeatedMirrorOutput=1;
				end
			case {'lineWidth', 'lineWidthSlider'}		% position 3 in the config
				state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 3)=state.blaster.lineWidth;
				if needUpdate
					state.internal.needNewPcellRepeatedOutput=1;
					state.internal.needNewRepeatedMirrorOutput=1;
				end
			case {'linePower1', 'linePowerSlider1'}		% position 4 in the config
				state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 4)=state.blaster.linePower1;
				if needUpdate
					state.internal.needNewPcellRepeatedOutput=1;
				end
			case {'linePower2', 'linePowerSlider2'}		% position 5 in the config
				state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 5)=state.blaster.linePower2;
                try                    
                    disp(['Laser 2 blaster power = ' ...
                            num2str(state.pcell.pcellPowerCal2(1+round(state.blaster.linePower2)))...
                            ' mW   '...
                           num2str(100*state.pcell.pcellPowerCal2(1+round(state.blaster.linePower2))/max(state.pcell.pcellPowerCal2))...
                           '% max']);
                    state.pcell.pcellBlasterPower2=state.pcell.pcellPowerCal2(1+round(state.blaster.linePower2));
                    updateHeaderString('state.pcell.pcellBlasterPower2');
                catch
                    state.pcell.pcellBlasterPower2=-1;
                    updateHeaderString('state.pcell.pcellBlasterPower2');
                end

				if needUpdate
					state.internal.needNewPcellRepeatedOutput=1;
				end
			case 'lineTilerActive'		% position 6 in the config
				state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, 6)=state.blaster.lineTilerActive;
				if needUpdate
					state.internal.needNewPcellRepeatedOutput=1;
				end
		end
		if needUpdate
			applyChangesToOutput;
		end
	catch
		disp(['genericConfig_Callback: ' lasterr]);
	end
	set(hObject, 'Enable', 'on');	


function varargout = genericUpdate_Callback(hObject, eventdata, handles, varargin)
	set(hObject, 'Enable', 'off');	
	try
		genericCallback(hObject);
		handleChange;
	catch
		disp(['blaster: ' lasterr]);
	end
	set(hObject, 'Enable', 'on');	

function varargout = commonParam_Callback(hObject, eventdata, handles, varargin)
	set(hObject, 'Enable', 'off');	
	try
		genericCallback(hObject);
		handleChange;
	catch
		disp(['blaster: ' lasterr]);
	end
	set(hObject, 'Enable', 'on');	
	

function displayPos_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	global state

	if state.blaster.displayPos>length(state.blaster.indexList)
		state.blaster.displayPos=length(state.blaster.indexList)+1;
		updateGuiByGlobal('state.blaster.displayPos');
		state.blaster.XList(state.blaster.displayPos)=0;
		state.blaster.YList(state.blaster.displayPos)=0;
	end
	updatePositionDisplay;

function pickPosition_Callback(hObject, eventdata, handles)
	global state
	if state.internal.status==2 | state.internal.status==3
		return
	end
	set(hObject, 'Enable', 'off');
	try
		setStatusString('Click on blast location');
		if state.internal.selectionChannel<=3
			figure(state.internal.GraphFigure(state.internal.selectionChannel));
		elseif state.internal.selectionChannel==4
			figure(state.internal.refFigure);
		end
		[x,y]=ginput(1);


       	c = cos(state.acq.scanRotation*pi/180);
    	s = sin(state.acq.scanRotation*pi/180);


%         offsetSmallMirror=(1-0.7*abs(s))*0.025;        offset in fractions of 1 line to center the blast when zoom=35
%         offsetBigMirror=(abs(s))*0.04;            offset in fractions of 1 line to center the blast when zoom=35

% 		offsetSmallMirror=0.105;        offset in fractions of 1 line to center the blast when zoom=35
%         offsetBigMirror=0.11;            offset in fractions of 1 line to center the blast when zoom=35
% 
%     	offsetX = abs(c*offsetSmallMirror + s*offsetBigMirror);  at 0, all small at 90 all big
%  		offsetY = abs(s*offsetSmallMirror + c*offsetBigMirror);  at 0, all big, at 90 all small
		
%  		offsetY = abs(offsetBigMirror * c - s * offsetSmallMirror);
   
%         offsetX=.11;
		offsetX=0;
   		offsetY=0;
   
		state.blaster.indexXList(state.blaster.displayPos)=round(x);
		state.blaster.indexYList(state.blaster.displayPos)=round(y);

        x=round(x - offsetX*state.acq.pixelsPerLine);
		y=round(y - offsetY*state.acq.linesPerFrame);

%        [offsetSmallMirror offsetBigMirror]
%        [offsetX offsetY]
%        [x y]
%        disp('blha')
  
		timeFromStartOfLine = state.acq.lineDelay+state.acq.mirrorLag + 1000*state.acq.pixelTime*(x-1);

        index = 1 + round(y-1) * state.internal.lengthOfXData ...
            + round(state.internal.lengthOfXData*timeFromStartOfLine/(1000*state.acq.msPerLine));
        
        
   %    state.internal.fractionStart + x*state.internal.fractionPerPixel));
		
      %  [state.acq.rotatedMirrorData(index,1)        state.acq.rotatedMirrorData(index,2)]
        
		state.blaster.X = state.acq.rotatedMirrorData(index,1);
		state.blaster.Y = state.acq.rotatedMirrorData(index,2);
		
		updateGuiByGlobal('state.blaster.Y');
		updateGuiByGlobal('state.blaster.X');
		
		state.blaster.XList(state.blaster.displayPos)=state.blaster.X;
		state.blaster.YList(state.blaster.displayPos)=state.blaster.Y;
		state.blaster.indexList(state.blaster.displayPos)=index;
		
		setStatusString('');

		if state.blaster.active & any(state.blaster.displayPos==state.blaster.allConfigs{state.blaster.currentConfig, 2}(:,1))	% the position is used in the current config
			state.internal.needNewRepeatedMirrorOutput=1;
			applyChangesToOutput;
		end
	catch
		setStatusString('ERROR');
		disp(['blaster selection : ' lasterr]);
	end
	
	try
		updateReferenceImage
	catch
		disp(['blaster selection (updateReferenceImage): ' lasterr]);
	end
	set(hObject, 'Enable', 'on');
	
function X_Callback(hObject, eventdata, handles)

function Y_Callback(hObject, eventdata, handles)

function lineWidth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function lineWidth_Callback(hObject, eventdata, handles)

function linePower1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function linePower1_Callback(hObject, eventdata, handles)

function config_CreateFcn(hObject, eventdata, handles)
	if ispc
        set(hObject,'BackgroundColor','white');
	else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
	end

function config_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	updateBlasterConfigDisplay
	
function line_CreateFcn(hObject, eventdata, handles)
	if ispc
        set(hObject,'BackgroundColor','white');
	else
        set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
	end

function line_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	global state
	if state.blaster.line > size(state.blaster.allConfigs{state.blaster.displayConfig, 2}, 1)
		state.blaster.line=size(state.blaster.allConfigs{state.blaster.displayConfig, 2}, 1)+1;
		updateGuiByGlobal('state.blaster.line');
		state.blaster.allConfigs{state.blaster.displayConfig, 2}(end+1, :)=	[state.blaster.linePos 0 state.blaster.lineWidth state.blaster.linePower1 state.blaster.linePower2 0];
		updateBlasterConfigLineDisplay;
		state.blaster.maxLine=state.blaster.line;
		updateGuiByGlobal('state.blaster.maxLine')
	else
		updateBlasterConfigLineDisplay;
	end

function deleteLine_Callback(hObject, eventdata, handles)
	global state
	if size(state.blaster.allConfigs{state.blaster.displayConfig, 2},1) > 1
		state.blaster.allConfigs{state.blaster.displayConfig, 2}(state.blaster.line, :)=[];
		state.blaster.maxLine=size(state.blaster.allConfigs{state.blaster.displayConfig, 2},1);
		updateGuiByGlobal('state.blaster.maxLine')
		if state.blaster.line > size(state.blaster.allConfigs{state.blaster.displayConfig, 2}, 1)
			state.blaster.line=size(state.blaster.allConfigs{state.blaster.displayConfig, 2}, 1);
			updateGuiByGlobal('state.blaster.line');
		end
		updateBlasterConfigLineDisplay;
	else	
		beep;
		setStatusString('CAN''T DELETE ONLY LINE');
	end
		
function maxLine_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
function maxLine_Callback(hObject, eventdata, handles)


		
function configName_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	global state
	state.blaster.allConfigs{state.blaster.displayConfig,1}=state.blaster.displayConfigName;
	makeBlasterConfigMenu

function updatePositionDisplay
	global state
	state.blaster.X = state.blaster.XList(state.blaster.displayPos);
	updateGuiByGlobal('state.blaster.X');
	state.blaster.Y = state.blaster.YList(state.blaster.displayPos);
	updateGuiByGlobal('state.blaster.Y');

function handleChange
	global state 	

	if state.blaster.active
		state.internal.needNewPcellRepeatedOutput=1;
		state.internal.needNewRepeatedMirrorOutput=1;
		applyChangesToOutput;
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS RELATED TO THE TILER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectArea_Callback(hObject, eventdata, handles)
	global state
	if state.internal.status==2 | state.internal.status==3
		return
	end
	set(hObject, 'Enable', 'off');
	try
		figure(state.internal.selectionChannel);
		k = waitforbuttonpress;
		setStatusString('Select Region');
	
		if isempty(findobj(gcf, 'Type', 'axes'))
			disp('*** NO axes***');
			return
		end
				
		point1 = get(gca,'CurrentPoint');    % button down detected
		finalRect = rbbox;                   % return figure units
		setStatusString('');
	
		point2 = get(gca,'CurrentPoint');    % button up detected
		point1 = point1(1,1:2);              % extract x and y
		point2 = point2(1,1:2);
		p1 = min(point1,point2);             % calculate locations
		offset = abs(point1-point2);         % and dimensions
		x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
		y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
		x=round(x);
		y=round(y);
		state.blaster.tilerAreaX0=x(1);
		state.blaster.tilerAreaX1=x(2);
		state.blaster.tilerAreaY0=y(1);
		state.blaster.tilerAreaY1=y(3);
		remakeTiler;
	catch
		setStatusString('ERROR');
		disp(['tiler selection : ' lasterr]);
	end
	set(hObject, 'Enable', 'on');
		
		
function tileCounter_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
function tileCounter_Callback(hObject, eventdata, handles)
function tileX_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
function tileX_Callback(hObject, eventdata, handles)
function tileY_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
function tileY_Callback(hObject, eventdata, handles)

function remakeTiler_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	remakeTiler;

function resetTiler_Callback(hObject, eventdata, handles)
	remakeTiler;
	
function reviewTile_Callback(hObject, eventdata, handles)
	genericCallback(hObject);
	global state

	state.blaster.tileX = state.acq.rotatedMirrorData(state.blaster.tilerIndexList(state.blaster.reviewTile), 1);
	state.blaster.tileY = state.acq.rotatedMirrorData(state.blaster.tilerIndexList(state.blaster.reviewTile), 2);
	updateGUIByGlobal('state.tiler.tileX');
	updateGUIByGlobal('state.tiler.tileY');

	if state.blaster.posUsingTiler
		state.blaster.indexList(state.blaster.posUsingTiler) = state.blaster.tilerIndexList(state.blaster.reviewTile);
		state.blaster.XList(state.blaster.posUsingTiler) = state.blaster.tileX;
		state.blaster.YList(state.blaster.posUsingTiler) = state.blaster.tileY;

		if state.blaster.displayPos == state.blaster.posUsingTiler
			state.blaster.X=state.blaster.tileX;
			state.blaster.Y=state.blaster.tileY;
			updateGUIByGlobal('state.blaster.X');
			updateGUIByGlobal('state.blaster.Y');
		end
		updateReferenceImage;
	end

	if length(state.tiler.acqList)>=state.tiler.reviewTile
		try
			reviewPhysAcq(state.tiler.acqList(state.tiler.reviewTile));
			reviewFluorData(state.tiler.acqList(state.tiler.reviewTile));
		end
	end


% --- Executes during object creation, after setting all properties.
function linePower2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linePower2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function linePower2_Callback(hObject, eventdata, handles)
% hObject    handle to linePower2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of linePower2 as text
%        str2double(get(hObject,'String')) returns contents of linePower2 as a double


% --- Executes during object creation, after setting all properties.
function linePowerSlider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to linePowerSlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function linePowerSlider2_Callback(hObject, eventdata, handles)
% hObject    handle to linePowerSlider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


