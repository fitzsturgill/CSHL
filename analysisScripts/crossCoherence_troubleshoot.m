    Fs = 610;
    dT = 1/Fs;
    t = linspace(0, 119.9984, 73200);
    freq = 2;
    stimTime = 40;
    responseDuration = 10;
    responseAmp = 3;
    sinAmp = 2;
    baseline = 4;
    y = sin(2 * pi * freq * t)*sinAmp/2 + sinAmp + baseline;
    response = linspace(1, responseAmp, responseDuration/2 *  Fs);
    response = [response fliplr(response)]; % make a triangle
    envelope = ones(size(t));
    startPoint = bpX2pnt(stimTime, Fs, 0);
    envelope(startPoint:startPoint + length(response) - 1) = response;
    t = t(:);
    envelope = envelope(:);
    y = y(:);
    y = y .* envelope;
    trials = 20;
    data1 = repmat(y, 1, trials);
    data2 = data1;
    ensureFigure('test', 1); plot(t, y);
    %%
    noiseAmp = 1;
    nSamples = length(data1);
    for counter = 1:trials
        noise1 = pinknoise(nSamples) * noiseAmp;
        noise2 = pinknoise(nSamples) * noiseAmp;
        data1(:, counter) = y + noise1(:);
        data2(:, counter) = y + noise2(:);        
    end
        
    
    
    %%
    params.Fs = Fs;
    params.trialave = 1;
    params.err = [2 0.05];
    params.tapers = [3 5];
    params.pad = 1;
    params.fpass = [0 20];
    
    movingwin = [5 1];
    
        cxcg = struct(...
        'C', [],...
        'phi', [],...
        'S12', [],...
        'S1', [],...
        'S2', [],...
        't', [],...    
        'f', [],...
        'Fs', Fs...
        );


    [cxcg.C, cxcg.phi, cxcg.S12, cxcg.S1, cxcg.S2, cxcg.t, cxcg.f] = cohgramc(data1, data2, movingwin, params);
    
    %% plot the data
    ensureFigure('cc_troubleshoot', 1);
    subplot(2,2,1); plot(t, y); title('before noise'); set(gca, 'XLim', [0 120]);
    subplot(2,2,2); plot(t, data1(:,1), 'g'); hold on; plot(t, data2(:, 1), 'r');
    title('after noise'); set(gca, 'XLim', [0 120]);
    subplot(2,2,3); 
    
    image(cxcg.t, cxcg.f, cxcg.C, 'CDataMapping', 'Scaled');
    colormap('jet');
    set(gca, 'Clim', [min(min(cxcg.C)), max(max(cxcg.C))]);
    title('cross coherence');