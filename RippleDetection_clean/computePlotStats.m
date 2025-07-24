function plotStatsFinal =  computePlotStats(events, fs)



for ch = 1:size(events.cov,2)
    
    
    nRip = length(events.cov{ch});
    centered = ones(nRip,1);
    
    for rip = 1:nRip

        %Ripple Amplitude
        plotStats.rippleAmp{ch}(rip) = events.rippleAmp{ch}(rip);

        %Ripple zero cross
        center = round(size(events.raw{ch}(rip, :),2)/2);
        rippleBand = events.rippleband{ch}(rip, :);
        RBzscore = events.RBzscore{ch}(rip,:);

        score = Inf;
        winStart = 0;
        while score > 0.75 && winStart < 200
            score = abs(RBzscore(center-winStart));
            winStart = winStart + 1;
        end

        score = Inf;
        winEnd = 0;
        while score > 0.5 && winEnd < 2000 
            score = abs(RBzscore(center+winEnd));
            winEnd = winEnd + 1;
        end

        plotStats.duration{ch}(rip) = (winStart + winEnd) * (1/fs) * 1000; %ms
        plotStats.window{ch}(rip,:) = [center-winStart, center+winEnd];

        window = center-winStart:center+winEnd;
        winData = rippleBand(window);
        zcs = find(winData(1:end-1).*winData(2:end)<0); %zero crossings

        freq = round((length(zcs)/2) / (length(window)/fs)*10)/10;

        plotStats.oscFreq{ch}(rip) = freq;

        if plotStats.duration{ch}(rip) < 3
            centered(rip) = 0;
        end


    end
    
    centered(centered==0)=1;
    centered = logical(centered);
    if length(centered) > 1
    %%%%% Patch - need to sort out ripple centering issue
        plotStatsFinal.InterRipPeriod{ch} = diff(events.goodRipples{ch}(centered));
        plotStatsFinal.duration{ch} = plotStats.duration{ch}(centered);
        plotStatsFinal.oscFreq{ch} = plotStats.oscFreq{ch}(centered);
        plotStatsFinal.rippleAmp{ch} = plotStats.rippleAmp{ch}(centered);
        plotStatsFinal.centeredInd{ch} = centered;
        plotStatsFinal.window{ch} = plotStats.window{ch}(centered,:);
    end    


end


return