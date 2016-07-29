function mcFigureKeyPressFunction
	global state
	val = double(get(gcbo,'CurrentCharacter'));

	if ~isnumeric(val) || isempty(val);
		return
	end
	try
		switch val
            
%             should have written mcZoom to accept a 2 element vector for the factors for shifting/zooming for  cases
            case 97 % a- xAxis Auto
                mcZoom('X', 'auto');
            case 115 % s- X Axis Standard [-.02 0.02]
                mcZoom('X', 'standard');
            case 100 % d- X zoom in
                mcZoom('X', 'in');
            case 102 % f- X zoom out
                mcZoom('X', 'out');
            case 113 % q- Y auto
                mcZoom('Y', 'auto');                
            case 119 % w - Y standard
                mcZoom('Y', 'standard');                
            case 101 % e - Y in
                mcZoom('Y', 'in');                
            case 114 % r - Y out
                mcZoom('Y', 'out');
            case 29 % right- shift x axes 200 to right
                mcZoom('X', 'right');
            case 28 % left - shift left
                mcZoom('X', 'left');
            case 30 % up- shift up by .005
                mcZoom('Y', 'up');
            case 31 % down- shift down
                mcZoom('Y', 'down');
%             case 109 % m-   advance fileCounter
%                 state.mcViewer.fileCounter = min(state.mcViewer.fileCounter +1, state.mcViewer.tsNumberOfFiles);
%                 updateGUIByGlobal('state.mcViewer.fileCounter');
%                 mcFlipTimeSeries;
%             case 110 %n - step back fileCounter
%                 state.mcViewer.fileCounter = max(state.mcViewer.fileCounter - 1, 1);
%                 updateGUIByGlobal('state.mcViewer.fileCounter');
%                 mcFlipTimeSeries;
            case 83 %S - auto for mcChannels, standard for aux channels or mcChannels that are differentially filtered via channel control GUIs (y axis only)
                mcZoom('Y', 'special');
		otherwise
		end
	catch
 		lasterr
    end

