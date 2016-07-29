function scaleDataViewerY(handles, xFactor, yFactor)
%
%
%   Created: 4/5/10 - SRO
%   Modified:

% Get handles for DataViewer plot objects
temp = getappdata(handles.hDataViewer,'handlesPlot');
hAllAxes = temp(1,:);
hPlotLines = temp(2,:);
hRasters = temp(3,:);
dvplot = getappdata(handles.hDataViewer,'dvplot');
PlotVectorOn = dvplot.pvOn;

% Display data
globalLim = [-0.001 0.001];
for i = PlotVectorOn'  % Must be row vector for this notation to work
    yAxes = get(hPlotLines(i),'YLim');
    
    
    deltaY = yAxes(2) - yAxes(1);
    deltaY2 = deltaY * yFactor;
    yMaxRatio = 
    
    
    
    if 0
        yLim = [min([temp globalLim(1)]) max([temp globalLim(2)])];
    else
        yLim = [min(temp) max(temp)];
    end
    set(hAllAxes(i),'YLim',yLim);
    % Update ticks
    setAxisTicks(hAllAxes(i));
end