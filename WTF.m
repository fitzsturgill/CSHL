% %     spacing = 0.1; % 0.1s steps to determine optimum window
%     assert(Fs == 20, 'Sample rate assumed to be 20');
%     commonCueWindow = [0 min(csWindow(:,2))]; % in case you've changed the delay across sessions, base windows upon smallest delay
% %     find peak us response
%     avgData = phAverageFromTE(TE, csPlusTrials & hitTrials, 1, 'FluorDataField', 'ZS', 'zeroTimes', TE.Cue, 'window', commonCueWindow); %high value, reward
% %     ensureFigure('test', 1);
% %     plot(avgData.xData, avgData.Avg);
%     [~, mix] = max(avgData.Avg);
% %     mix = fix(mix/2)*2; % make even
%     prePoints = (1:2:mix)';
%     postPoints = mix+1:2:bpX2pnt(commonCueWindow(2), 20, 0);
%     % make matrix of pre, post points, and D' for each pairing
%     preMatrix = repmat(prePoints, 1, length(postPoints));
%     postMatrix = repmat(postPoints, length(prePoints), 1);
%     dMatrix = zeros(length(prePoints), length(postPoints), 2);    
%     
%     % collect phWindow for commonCuewindow for each channel
%     for ch = 1:2        
%         [phData_plus, ~] = phAlignedWindow(TE, csPlusTrials & hitTrials, ch, 'zeroTimes', TE.Cue, 'window', commonCueWindow, 'FluorDataField', 'ZS');
%         [phData_minus, ~] = phAlignedWindow(TE, csMinusTrials & CRTrials, ch, 'zeroTimes', TE.Cue, 'window', commonCueWindow, 'FluorDataField', 'ZS');
%         for counter = 1:numel(preMatrix)
%             [D, P] = rocarea(nanmean(phData_plus(:,preMatrix(counter):postMatrix(counter)), 2), nanmean(phData_minus(:,preMatrix(counter):postMatrix(counter)), 2), 'scale');
%             dMatrix(counter + (ch - 1) * numel(preMatrix)) = D;
%         end
%     end