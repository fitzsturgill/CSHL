function out = mcFilter(data, lowCutoff, highCutoff, Fs)
% function out = mcFilter(data, lowCutoff, highCutoff, Fs)
% low level filter function
% data- column vector or array with channels oriented in columns
% lowCutoff- cutoff freq. for lowpass filter
% highCutoff- cutoff freq. for highpass filter
% Fs- sampling frequency

    global state
    if nargin <4
        Fs = state.mcViewer.sampleRate;
    end

    % Create butterworth filter
    Fs = Fs * 1000;
      
    lowCutoff = lowCutoff*2/Fs;
    highCutoff = highCutoff*2/Fs;

    % note-   5 pole Bessel filter in Matlab used in Frohlich and McCormick  
    if lowCutoff && highCutoff
        [b, a] = butter(5, highCutoff, 'high');
        out = filtfilt(b,a,data);
        [b, a] = butter(5, lowCutoff, 'low');
        out = filtfilt(b,a,out);        
    elseif highCutoff
        [b, a] = butter(5, highCutoff, 'high');
        out = filtfilt(b,a,data);        
    elseif lowCutoff
        [b, a] = butter(5, lowCutoff, 'low');         
        out = filtfilt(b,a,data);
    else
        out=data;
    end
    
    

    
    
    
    
    
% function out = mcFilter(data, lowCutoff, highCutoff, Fs)
% low level filter function
% data- column vector or array with channels oriented in columns
% lowCutoff- cutoff freq. for lowpass filter
% highCutoff- cutoff freq. for highpass filter
% Fs- sampling frequency

%     global state
%     if nargin <4
%         Fs = state.mcViewer.sampleRate;
%     end
% 
%     % Create butterworth filter
%     Fs = Fs * 1000;
%     
%     lowStopBand = min((lowCutoff + 20) * 2/Fs, 0.9999);
%     highStopBand = max((highCutoff - 20) * 2/Fs, 0.0001);
% %     lowStopBand = lowCutoff * (3/4) * 2/Fs;
% %     highStopBand = highCutoff * (5/4) * 2/Fs;    
%     lowCutoff = lowCutoff*2/Fs;
%     highCutoff = highCutoff*2/Fs;
%     if lowCutoff && highCutoff
%         Wn = [highCutoff lowCutoff];
%         Wp = [highCutoff lowCutoff];
%         Ws = [highStopBand lowStopBand];
%     elseif lowCutoff
%         Wn = lowCutoff;
%         Wp = lowCutoff;
%         Ws = lowStopBand;
%     elseif highCutoff
%         Wn = highCutoff;
%         Wp = highCutoff;
%         Ws = highStopBand;
%     else
%         out = data;
%         return
%     end
% 
%     [n, Wn] = buttord(Wp, Ws, 3, 20); 
%     if lowCutoff && highCutoff
%         [b, a] = butter(n, Wn);        
%     elseif highCutoff
%         [b, a] = butter(5, highCutoff, 'high');
%     else
%         [b, a] = butter(5, lowCutoff, 'low');   % note-   5 pole Bessel filter in Matlab used in Frohlich and McCormick        
%     end
%     
%     % filter data
%     out = filtfilt(b,a,data);    