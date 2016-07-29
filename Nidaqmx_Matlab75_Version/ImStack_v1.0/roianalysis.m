function varargout = roianalysis(varargin)
% ROIANALYSIS Application M-file for roianalysis.fig
%    FIG = ROIANALYSIS launch roianalysis GUI.
%    ROIANALYSIS('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 15-Jul-2002 15:25:31

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
function varargout = roiPopupmenu_Callback(h, eventdata, handles, varargin)
global state gh
value=get(h,'Value');
state.imageProc.roiAnalysis.roiFrame = state.imageProc.cell.currentFrame{value};
state.imageProc.roiAnalysis.roiEnd = size(state.imageProc.cell.currentImage{value},3);
state.imageProc.roiAnalysis.minPixel= state.imageProc.cell.lowPixelValue{value};
state.imageProc.roiAnalysis.maxPixel=state.imageProc.cell.highPixelValue{value};
state.imageProc.roiAnalysis.imageData = state.imageProc.cell.currentImage{value}(:,:,state.imageProc.roiAnalysis.roiFrame);
set(gh.roiAnalysis.roiAxes,'CLim',[state.imageProc.roiAnalysis.minPixel state.imageProc.roiAnalysis.maxPixel],...
    'XLim',[1 size(state.imageProc.roiAnalysis.imageData,2)],'YLim',[1 size(state.imageProc.roiAnalysis.imageData,1)]);

updateGUIByGlobal('state.imageProc.roiAnalysis.roiFrame');
updateGUIByGlobal('state.imageProc.roiAnalysis.roiEnd');
updateGUIByGlobal('state.imageProc.roiAnalysis.minPixel');
updateGUIByGlobal('state.imageProc.roiAnalysis.maxPixel');
if ishandle(state.imageProc.roiAnalysis.imageHandle)
    set(state.imageProc.roiAnalysis.imageHandle,'CData',state.imageProc.roiAnalysis.imageData);
    set(gh.roiAnalysis.figure1,'ColorMap', makeColorMap('gray',8));
else
    state.imageProc.roiAnalysis.imageHandle=image('CData',state.imageProc.roiAnalysis.imageData,...
        'CDataMapping','Scaled','Parent',gh.roiAnalysis.roiAxes);
    set(gh.roiAnalysis.figure1,'ColorMap', makeColorMap('gray',8));
end



% --------------------------------------------------------------------
function varargout = selectROI_Callback(h, eventdata, handles, varargin)
global gh state
rectCoords=getRect(gh.roiAnalysis.roiAxes);
patchCoords=convertRectToPatch(rectCoords);
if ishandle(state.imageProc.roiAnalysis.patchHandle)
    delete(state.imageProc.roiAnalysis.patchHandle);
	state.imageProc.roiAnalysis.roiAngle=0;
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiAngle');
end
state.imageProc.roiAnalysis.patchHandle=patch(patchCoords(1,:),patchCoords(2,:),[1 1 0], ...
    'AlphaDataMapping','none','FaceAlpha',.6,'Parent',gh.roiAnalysis.roiAxes,'ButtonDownFcn','roiPatchButtonDwn');
state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');

% Do analysis
roiStats;



% --------------------------------------------------------------------
function varargout = roiAngle_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
% Rotate object about its center....
if ishandle(state.imageProc.roiAnalysis.patchHandle)
    state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
    state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
    state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
    state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
    state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
    
	rotate(state.imageProc.roiAnalysis.patchHandle,[0 0 1],...
        state.imageProc.roiAnalysis.roiAngle-state.imageProc.roiAnalysis.lastAngle,state.imageProc.roiAnalysis.patchCenter );

    state.imageProc.roiAnalysis.lastAngle=state.imageProc.roiAnalysis.roiAngle;
	
    state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
    state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
    state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
    state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
    state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
end

% Do analysis
roiStats;


% --------------------------------------------------------------------
function varargout = roiAngleSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
% Rotate object about its center....
if ishandle(state.imageProc.roiAnalysis.patchHandle)
    state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
    state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
    state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
    state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
    state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
	
  rotate(state.imageProc.roiAnalysis.patchHandle,[0 0 1],...
        state.imageProc.roiAnalysis.roiAngle-state.imageProc.roiAnalysis.lastAngle,state.imageProc.roiAnalysis.patchCenter );
	state.imageProc.roiAnalysis.lastAngle=state.imageProc.roiAnalysis.roiAngle;
    
    state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
    state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
    state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
    state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
    state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
    updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
    updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');
end
% Do analysis
roiStats;





% --------------------------------------------------------------------
function varargout = roixpos_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
if ishandle(state.imageProc.roiAnalysis.patchHandle)
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
	CurrentroiXpos=mean(state.imageProc.roiAnalysis.XPatch);
	difference=state.imageProc.roiAnalysis.roixpos-CurrentroiXpos;
	set(state.imageProc.roiAnalysis.patchHandle,'XData',state.imageProc.roiAnalysis.XPatch+difference);
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
end

% Do analysis
roiStats;

% --------------------------------------------------------------------
function varargout = roixposSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
if ishandle(state.imageProc.roiAnalysis.patchHandle)
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
	CurrentroiXpos=mean(state.imageProc.roiAnalysis.XPatch);
	difference=state.imageProc.roiAnalysis.roixpos-CurrentroiXpos;
	set(state.imageProc.roiAnalysis.patchHandle,'XData',state.imageProc.roiAnalysis.XPatch+difference);
	state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
end

% Do analysis
roiStats;


% --------------------------------------------------------------------
function varargout = roiypos_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
if ishandle(state.imageProc.roiAnalysis.patchHandle)
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
	CurrentroiYpos=mean(state.imageProc.roiAnalysis.YPatch);
	difference=state.imageProc.roiAnalysis.roiypos-CurrentroiYpos;
	set(state.imageProc.roiAnalysis.patchHandle,'YData',state.imageProc.roiAnalysis.YPatch+difference);
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
end

% Do analysis
roiStats;




% --------------------------------------------------------------------
function varargout = roiyposSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
if ishandle(state.imageProc.roiAnalysis.patchHandle)
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
	CurrentroiYpos=mean(state.imageProc.roiAnalysis.YPatch);
	difference=state.imageProc.roiAnalysis.roiypos-CurrentroiYpos;
	set(state.imageProc.roiAnalysis.patchHandle,'YData',state.imageProc.roiAnalysis.YPatch+difference);
	state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
end

% Do analysis
roiStats;


% --------------------------------------------------------------------
function varargout = roiFrame_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global state gh
value=get(gh.roiAnalysis.roiPopupmenu,'Value');
try
    state.imageProc.roiAnalysis.imageData = state.imageProc.cell.currentImage{value}(:,:,state.imageProc.roiAnalysis.roiFrame);
    if ishandle(state.imageProc.roiAnalysis.imageHandle)
        set(state.imageProc.roiAnalysis.imageHandle,'CData',state.imageProc.roiAnalysis.imageData);
    else
        state.imageProc.roiAnalysis.imageHandle=image('CData',state.imageProc.roiAnalysis.imageData,...
            'CDataMapping','Scaled','Parent',gh.roiAnalysis.roiAxes);
    end
end

% Do analysis
roiStats;


% --------------------------------------------------------------------
function varargout = roiFrameSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global state gh
value=get(gh.roiAnalysis.roiPopupmenu,'Value');
try
    state.imageProc.roiAnalysis.imageData = state.imageProc.cell.currentImage{value}(:,:,state.imageProc.roiAnalysis.roiFrame);
    if ishandle(state.imageProc.roiAnalysis.imageHandle)
        set(state.imageProc.roiAnalysis.imageHandle,'CData',state.imageProc.roiAnalysis.imageData);
    else
        state.imageProc.roiAnalysis.imageHandle=image('CData',state.imageProc.roiAnalysis.imageData,...
            'CDataMapping','Scaled','Parent',gh.roiAnalysis.roiAxes);
    end
end
 
% Do analysis
roiStats;

% --------------------------------------------------------------------
function varargout = roiEnd_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiEndSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiSum_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiMean_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiMedian_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiMax_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiMin_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = maxPixel_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
set(gh.roiAnalysis.roiAxes,'CLim',[state.imageProc.roiAnalysis.minPixel state.imageProc.roiAnalysis.maxPixel]);



% --------------------------------------------------------------------
function varargout = maxPixelSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
set(gh.roiAnalysis.roiAxes,'CLim',[state.imageProc.roiAnalysis.minPixel state.imageProc.roiAnalysis.maxPixel]);



% --------------------------------------------------------------------
function varargout = minPixel_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
set(gh.roiAnalysis.roiAxes,'CLim',[state.imageProc.roiAnalysis.minPixel state.imageProc.roiAnalysis.maxPixel]);



% --------------------------------------------------------------------
function varargout = minPixelSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);
global gh state
set(gh.roiAnalysis.roiAxes,'CLim',[state.imageProc.roiAnalysis.minPixel state.imageProc.roiAnalysis.maxPixel]);



% --------------------------------------------------------------------
function varargout = loadROI_Callback(h, eventdata, handles, varargin)
global gh state
disp('Data Written to array called ROIStats.');
% disp('ROIStats = [Sum Mean Median Min Max].');
stats=[state.imageProc.roiAnalysis.roiSum state.imageProc.roiAnalysis.roiMean ...	
    state.imageProc.roiAnalysis.roiStd state.imageProc.roiAnalysis.roiMedian ...
    state.imageProc.roiAnalysis.roiMax state.imageProc.roiAnalysis.roiMin ...
	state.imageProc.roiAnalysis.roiNumber];
	assignin('base','ROIStats',stats);
try
	eval(['ddepoke(state.imageProc.excellink.excelChannel, ' '''' ['r' num2str(state.imageProc.spine.row+1) 'c' num2str(1) ':r' num2str(state.imageProc.spine.row+1) 'c' num2str(7)] '''' ',stats'''' );']);
	state.imageProc.spine.row = state.imageProc.spine.row+1;
	updateGUIByGlobal('state.imageProc.spine.row');

end




% --------------------------------------------------------------------
function varargout = roiHistogram_Callback(h, eventdata, handles, varargin)
global gh state

if isempty(state.imageProc.roiAnalysis.roiData)
    beep;
    error('ROI Data not loaded');
end
[n,xout]=hist(state.imageProc.roiAnalysis.roiData,state.imageProc.roiAnalysis.roiBins);
figure('DoubleBuffer','On','NumberTitle','off','Name','ROI Pixel Histogram',...
    'Color','white','pos',[106   507   868   389]);
colormap(gray);
bar(xout,n);



% --------------------------------------------------------------------
function varargout = loadImage_Callback(h, eventdata, handles, varargin)
global state gh
value=get(gh.roiAnalysis.roiPopupmenu,'Value');

state.imageProc.roiAnalysis.roiFrame = state.imageProc.cell.currentFrame{value};
state.imageProc.roiAnalysis.roiEnd = size(state.imageProc.cell.currentImage{value},3);
state.imageProc.roiAnalysis.minPixel= state.imageProc.cell.lowPixelValue{value};
state.imageProc.roiAnalysis.maxPixel=state.imageProc.cell.highPixelValue{value};
state.imageProc.roiAnalysis.imageData = state.imageProc.cell.currentImage{value}(:,:,state.imageProc.roiAnalysis.roiFrame);
set(gh.roiAnalysis.roiAxes,'CLim',[state.imageProc.roiAnalysis.minPixel state.imageProc.roiAnalysis.maxPixel],...
    'XLim',[1 size(state.imageProc.roiAnalysis.imageData,2)],'YLim',[1 size(state.imageProc.roiAnalysis.imageData,1)]);

updateGUIByGlobal('state.imageProc.roiAnalysis.roiFrame');
updateGUIByGlobal('state.imageProc.roiAnalysis.roiEnd');
updateGUIByGlobal('state.imageProc.roiAnalysis.minPixel');
updateGUIByGlobal('state.imageProc.roiAnalysis.maxPixel');
if ishandle(state.imageProc.roiAnalysis.imageHandle)
    set(state.imageProc.roiAnalysis.imageHandle,'CData',state.imageProc.roiAnalysis.imageData);
    set(gh.roiAnalysis.figure1,'ColorMap', makeColorMap('gray',8));
else
    state.imageProc.roiAnalysis.imageHandle=image('CData',state.imageProc.roiAnalysis.imageData,...
        'CDataMapping','Scaled','Parent',gh.roiAnalysis.roiAxes);
    set(gh.roiAnalysis.figure1,'ColorMap', makeColorMap('gray',8));
end







% --------------------------------------------------------------------
function varargout = roiBins_Callback(h, eventdata, handles, varargin)
genericCallback(h)



% --------------------------------------------------------------------
function varargout = roiBinsSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = zoomOut_Callback(h, eventdata, handles, varargin)
global gh state
if ~isempty(state.imageProc.roiAnalysis.imageData)
	set(gh.roiAnalysis.roiAxes,'XLim',[1 size(state.imageProc.roiAnalysis.imageData,2)],...
		'YLim',[1 size(state.imageProc.roiAnalysis.imageData,1)]);
end




% --------------------------------------------------------------------
function varargout = zoomIn_Callback(h, eventdata, handles, varargin)
global gh
zoomIn(gh.roiAnalysis.roiAxes);



% --------------------------------------------------------------------
function varargout = roiNumber_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiWidth_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = roiWidthSlider_Callback(h, eventdata, handles, varargin)
genericCallback(h);



% --------------------------------------------------------------------
function varargout = drawROI_Callback(h, eventdata, handles, varargin)
global gh state
rectCoords=getRect(gh.roiAnalysis.roiAxes);
rectCoords(3:4) = [state.imageProc.roiAnalysis.roiWidth state.imageProc.roiAnalysis.roiWidth];
patchCoords=convertRectToPatch(rectCoords);
if ishandle(state.imageProc.roiAnalysis.patchHandle)
    delete(state.imageProc.roiAnalysis.patchHandle);
	state.imageProc.roiAnalysis.roiAngle=0;
	updateGUIByGlobal('state.imageProc.roiAnalysis.roiAngle');
end
state.imageProc.roiAnalysis.patchHandle=patch(patchCoords(1,:),patchCoords(2,:),[1 1 0], ...
    'AlphaDataMapping','none','FaceAlpha',.6,'Parent',gh.roiAnalysis.roiAxes,'ButtonDownFcn','roiPatchButtonDwn');
state.imageProc.roiAnalysis.XPatch=get(state.imageProc.roiAnalysis.patchHandle,'XData');
state.imageProc.roiAnalysis.YPatch=get(state.imageProc.roiAnalysis.patchHandle,'YData');
state.imageProc.roiAnalysis.patchCenter = [mean(state.imageProc.roiAnalysis.XPatch) mean(state.imageProc.roiAnalysis.YPatch) 0];
state.imageProc.roiAnalysis.roixpos=mean(state.imageProc.roiAnalysis.XPatch);
state.imageProc.roiAnalysis.roiypos=mean(state.imageProc.roiAnalysis.YPatch);
updateGUIByGlobal('state.imageProc.roiAnalysis.roiypos');
updateGUIByGlobal('state.imageProc.roiAnalysis.roixpos');

% Do analysis
roiStats;




% --------------------------------------------------------------------
function varargout = roiStd_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = timeSeriesAnalysis_Callback(h, eventdata, handles, varargin)
timeSeriesROI;



% --------------------------------------------------------------------
function varargout = plotTSeries_Callback(h, eventdata, handles, varargin)
if evalin('base','exist(''TimeSeriesMean'')')
	if evalin('base','iswave(TimeSeriesMean);')
		evalin('base','plot(TimeSeriesMean);');
		
	end
end
