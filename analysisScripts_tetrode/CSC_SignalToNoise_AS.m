        [fname, pname] = uiputfile('path', 'Choose CSC path...');
        if pname == 0            
            return
        else
            filepath = pname;
            disp(pname);
        end

rmsNoise = CSC_RMS('filepath', filepath);
ensureFigure('RMSNoise', 1);
plot(rmsNoise);
saveas(gcf, fullfile(filepath, 'RMSNoise.fig'));
saveas(gcf, fullfile(filepath, 'RMSNoise.jpg'));




