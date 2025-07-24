function [cleanedEvents, rejectVector] = rejectRipples(ripple, RejectParams, rippleWindow, ch, fs, IISflag, spikeMask, IIScheck)

smthKrnl    = RejectParams.bndParams.smthKrnl;
srchWinSz   = RejectParams.bndParams.srchWinSz;
scoreThresh = RejectParams.bndParams.scoreThresh;

if nargin < 7
    spikeMask = zeros(1,max(ripple.ind));
else
    spikeMask = spikeMask(:,ch);
end
    
center = size(ripple.events.raw,2)/2;
center = round(center);
rippleRange = center-rippleWindow : center+rippleWindow;
% iisRange = center-10*rippleWindow  : center+10*rippleWindow;

keepRipples = ones(1,length(ripple.ind));   
spikeOverlap = zeros(1,length(ripple.ind));   
rejectVector = zeros(9  ,length(ripple.ind));
for rip = 1:length(ripple.ind)
    
    rippleLoc = ripple.ind(rip);
    RBz       = max(abs(ripple.events.RBzscore(rip,center))); 
%     RBamp     = ripple.events.rippleAmp(rip);
%     RBHFratio = ripple.events.highRippleRatio(rip);
    HFz       = max(abs(ripple.events.HFzscore(rip,rippleRange)));


    if RejectParams.RBHFzscore(2)
        cond1 =  (RBz - HFz) < RejectParams.RBHFzscore(1);
    else
        cond1 = 0;
    end
        
    if RejectParams.RBzscore(2) == 1
        cond4 = RBz < RejectParams.RBzscore(1);
    elseif RejectParams.RBzscore(2) == 2
        cond4 =  (RBz - HFz) < RejectParams.RBHFzscore(1) && RBz < RejectParams.RBzscore(1);
    else 
        cond4 = 0;
    end
    
    if RejectParams.LFPthresh(2)
        cond5 =  any(abs(ripple.events.raw(rip,:)) > RejectParams.LFPthresh(1)); 
    else
        cond5 = 0;
    end
    
    if IISflag 
        startInd = rippleLoc-IIScheck;
        if startInd < 1; startInd = 1; end
        endInd = rippleLoc+IIScheck;
        if endInd > length(spikeMask); endInd = length(spikeMask); end

        spikeMaskCheck = spikeMask(startInd:endInd);
        if any(spikeMaskCheck)
            cond6 = 1;
            spikeOverlap(rip) = 1;
        else 
            cond6 = 0;
        end
        
    else
        cond6 = 0;
    end
    
    % calc duration 
%     score = Inf;
%     winStart = 0;
%     while score > 0.75 && winStart < fs
%         score = ripple.events.RBzscore(rip,center-winStart);
%         winStart = winStart + 1;
%     end
% 
%     score = Inf;
%     winEnd = 0;
%     while score > 0.75 && winEnd < fs 
%         score = ripple.events.RBzscore(rip,center+winEnd);
%         winEnd = winEnd + 1;
%     end

     RBzscore = ripple.events.RBzscore(rip,:);
     RBzscoreSmooth = smoothdata(RBzscore,2,'gaussian',smthKrnl); %apply different smoothing for edge detection compared to ripple thresholding
        

    score = Inf;
    winStart = 0;
    while score > scoreThresh && winStart < fs*2%200
        score = RBzscoreSmooth(center-winStart);
        winStart = winStart + 1;
    end
    range=(-round(smthKrnl*(srchWinSz/2)):round(smthKrnl*(srchWinSz/2)))+center-winStart;
    range(range>center)=[];
    range(range<1)=[];
    %[~,idx] = min(RBzscore(range));
    if all(isnan(RBzscore(range)))
        winStart= (center-range(end)) - find(~isnan(RBzscore(range(end):center)),1,'first');
    elseif winStart~=fs*2
        idx=[];
        bump=0;
        while isempty(idx)
            idx=find(RBzscore(range)<scoreThresh+bump,1,'last');
            bump=bump+0.1;
        end
        winStart = winStart - (idx - round(smthKrnl*(srchWinSz/2)));
    end
    
    score = Inf;
    winEnd = 0;
    while score > scoreThresh && winEnd < fs*2 %200
        score = RBzscoreSmooth(center+winEnd);
        winEnd = winEnd + 1;
    end
    range=(-round(smthKrnl*(srchWinSz/2)):round(smthKrnl*(srchWinSz/2)))+center+winEnd;
    nTrim=sum(range<center);
    range(range<center)=[];
    range(range>length(RBzscore))=[];
    %[~,idx]=min(RBzscore(range));
    if all(isnan(RBzscore(range)))
        winEnd=find(~isnan(RBzscore(center:range(1))),1,'last');
    elseif winEnd~=fs*2
        idx=[];
        bump=0;
        while isempty(idx)
            idx=find(RBzscore(range)<scoreThresh+bump,1,'first');
            bump=bump+0.1;
        end
        winEnd = winEnd + (idx - (round(smthKrnl*(srchWinSz/2)) - nTrim));
    end

    if winStart >= size(ripple.events.raw,2)/2; winStart = floor(size(ripple.events.raw,2)/2) - 1; end
    if winEnd >= size(ripple.events.raw,2)/2; winEnd = floor(size(ripple.events.raw,2)/2) - 1; end
    

    if isempty(winStart) || isempty(winEnd)
        duration = 0;
    else
        duration = (winStart + winEnd) * (1/fs) * 1000; %ms
    end

    if duration < RejectParams.minDuration(1)*fs || duration > fs%3*(1/100 * fs) %80 Hz is putative ripple frequency    
        cond7 = 1;
    else
        cond7 = 0;
    end
    
    if RejectParams.LFPdiff(2)
        
        scale = fs / 1000; %scale by sample rate
        
        LFPgrad = ripple.events.diffLFP(rip,:);
        LFPgrad = LFPgrad(center - winStart:center + winEnd);
        cond3 =  max(LFPgrad)*scale > RejectParams.LFPdiff(1);
    else
        cond3 = 0;
    end
    
    if RejectParams.sharpzscore(2)
        sharpZ = ripple.events.sharpzscore(rip,center - winStart:center + winEnd);
        cond2 =  any(sharpZ > RejectParams.sharpzscore(1));
    else
        cond2 = 0;
    end
    
    if RejectParams.LFPdiffzscore(2)
        LFPdiffzscore = ripple.events.diffLFPzscore(rip,center - winStart:center + winEnd);
        cond8 =  any(LFPdiffzscore > RejectParams.LFPdiffzscore(1));
    else
        cond8 = 0;
    end
    
    
    

    if cond1 || cond2 || cond3 || cond4 || cond5 || cond6 || cond7 || cond8
        keepRipples(rip) = 0;
    end
    
    
    % look at borderline rejected events.
    
    if cond1
        rejectVector(2,rip) = abs((RBz - HFz) - RejectParams.RBHFzscore(1));
    end
    if cond2
        rejectVector(3,rip) = abs(max(sharpZ) - RejectParams.sharpzscore(1));
    end
    if cond3
        rejectVector(4,rip) = abs(max(LFPgrad) - RejectParams.LFPdiff(1));
    end
    if cond4
        rejectVector(5,rip) = abs(RBz - RejectParams.RBzscore(1));
    end
    if cond5
        rejectVector(6,rip) = abs(max(abs(ripple.events.raw(rip,:))) - RejectParams.LFPthresh(1));
    end
    if cond6
        rejectVector(7,rip) = 1;
    end
    if cond7
        rejectVector(8,rip) = duration;
    end
    if cond8
        rejectVector(9,rip) = max(LFPdiffzscore) - RejectParams.LFPdiffzscore(1);
    end
    
end

% rejectVal = sum(rejectVector,1);
% pcntle = quantile(rejectVal,[0.01]);

% rejectVector(2:end,:) = rejectVector(1:8,:);
rejectVector(1,:) = ripple.ind;

keepRipples = logical(keepRipples);

% keepRipples = rejectVal < pcntle;   

fn = fieldnames(ripple.events);
for f = 1:numel(fn)
    fDat = ripple.events.(fn{f});
    DIM = min(size(fDat));
    if DIM > 1
        cleanedEvents.(fn{f}) = fDat(keepRipples,:);
    elseif DIM == 1
        cleanedEvents.(fn{f}) = fDat(keepRipples);
    elseif DIM == 0
        cleanedEvents.(fn{f}) = [];
    end
end



cleanedEvents.goodRipples = ripple.ind(keepRipples);

density = sum(keepRipples) / length(ripple.rippleband) * fs * 60;
% fprintf(['Kept ', num2str(sum(keepRipples) / length(ripple.ind)), ' of ripples on channel ', num2str(ch), '\n']) 
fprintf(['Kept ', num2str(sum(keepRipples)), ' ripples on channel ', num2str(ch), '\n']) 
fprintf(['Density ', num2str(sum(density)),'\n']) 
fprintf([ num2str(sum(spikeOverlap) / length(ripple.ind)), ' overlapped w spikes on channel ', num2str(ch), '\n\n']) 


    
            
        







return 