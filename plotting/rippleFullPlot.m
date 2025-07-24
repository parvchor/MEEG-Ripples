function rippleFullPlot(species, events, col, ch, rip_num, fs, rippleband)

plot_win = round(fs * 0.1); % +/- from events center, may want to zoom in more
plot_center = round(size(events.raw{ch}(rip_num, :),2)/2);
plot_range = plot_center-plot_win:plot_center+plot_win;
plot_win = round(((plot_range - plot_center) / fs) * 1000);
shift = col - 1;

% subtightplot(7,4, 1+shift, [0.04, 0.04])
% str = sprintf(['events Amp : %.1f (%2.0f %%ile)\n' ...
%                '> 200Hz Amp: %.1f (%2.0f %%ile)\n' ...
%                '25-50Hz Amp: %.1f (%2.0f %%ile)\n' ...
%                'RipHighRat : %.1f (%2.0f %%ile)\n' ...
%                'RipLowRat  : %.1f (%2.0f %%ile)\n' ...
%                'cov        : %.1f (%2.0f %%ile)\n'], ...
%                events.rippleAmp{ch}(rip_num),  100 * findPercentile(events.rippleAmp{ch},events.rippleAmp{ch}(rip_num)), ...
%                events.highFreqAmp{ch}(rip_num),100 * findPercentile(events.highFreqAmp{ch},events.highFreqAmp{ch}(rip_num)), ...
%                events.lowFreqAmp{ch}(rip_num), 100 * findPercentile(events.lowFreqAmp{ch},events.lowFreqAmp{ch}(rip_num)), ...
%                events.highRippleRatio{ch}(rip_num), 100 * findPercentile(events.highRippleRatio{ch},events.highRippleRatio{ch}(rip_num)), ...
%                events.lowRippleRatio{ch}(rip_num), 100 * findPercentile(events.lowRippleRatio{ch},events.lowRippleRatio{ch}(rip_num)), ...
%                events.cov{ch}(rip_num), 100 * findPercentile(events.cov{ch},events.cov{ch}(rip_num)));

% t = text(-100, 0, str, 'BackgroundColor', 'white');
% t.FontSize = 8;
% 
% xlim([-200 200])
% ylim([-10 10])
% axis off

%%
% subtightplot(7,4, 5+shift:4:12+shift, [0.04, 0.04])
subtightplot(7,4, 1+shift:4:8+shift, [0.04, 0.06])

yyaxis left
plot(plot_win, events.raw{ch}(rip_num, plot_range), 'k', 'LineWidth', 1.0);
hold on
ylim([min(events.raw{ch}(rip_num, plot_range)) max(events.raw{ch}(rip_num, plot_range))])
ylabel('LFP \muV')

yyaxis right
p2 = plot(plot_win, events.rippleband{ch}(rip_num,plot_range), 'g-', 'LineWidth', 2.5); hold on
p2.Color = [0 .7 0 0.5];    
% ylim([-10 10])    
% l = legend(p2, {sprintf('%i-%iHz amp', rippleband)});
% l.Location = 'northwest';
ylabel('rippleband \muV')
xlabel('time to ripple center [ms]')


% inds = find(events.goodRipples{ch} ~= events.goodRipples{ch}(rip_num));
% check = ((events.goodRipples{ch}(inds) > events.goodRipples{ch}(rip_num) - plot_win) & (events.goodRipples{ch}(inds) < events.goodRipples{ch}(rip_num) + plot_win));
% if any(check)
%     closeEvents = events.goodRipples{ch}(inds(check));
%     RelTime = closeEvents - events.goodRipples{ch}(rip_num);  
%     RelTimeRip = RelTime + plot_center;
%     p = plot(RelTime,events.raw{ch}(rip_num, RelTimeRip),'rs');
%     hold on;
%     p.MarkerFaceColor = 'r';
% end

vline(0)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0 .7 0];

%%
subtightplot(7,4, 9+shift:4:13+shift,   [0.04, 0.06])
% yyaxis right
% p3 = plot(-plot_win:plot_win, events.highFreqBand{ch}(rip_num,plot_range),  ...
%           'r-', 'LineWidth', 0.25, 'DisplayName', '200+ Hz'); hold on
% % ylim([-5 5])

% yyaxis left
% p4 = plot(-plot_win:plot_win, events.lowFreqBand{ch}(rip_num,plot_range), ...
%           'b-', 'LineWidth', 0.25, 'DisplayName', '25-50 Hz'); hold on
% % ylim([-10 10])
yyaxis right
ln1 = plot(plot_win, abs(hilbert(events.highFreqBand{ch}(rip_num,plot_range))), 'r-'); hold on;
% ln1.LineWidth = 2;
% ln2 = plot(-plot_win:plot_win, smoothdata(events.zscore100200{ch}(rip_num,plot_range), 'movmean',1), 'g-'); hold on;
switch species
    case 'human'
        ylabel('200+Hz amp \muV')
    case 'rodent'
        ylabel('300+Hz amp \muV')
end

yyaxis left
ln2 = plot(plot_win, abs(hilbert(events.band100200{ch}(rip_num,plot_range))), 'g-'); hold on;
ln2.LineWidth = 2;
ln2.Color = [0 0 0.6];

% legend([ln2, ln1], {'100-200 Hz','200+ Hz'})


vline(0)
ylabel('100-200Hz amp \muV')
% ylabel('z score')

ax = gca;
ax.YAxis(2).Color = ln1.Color; 
ax.YAxis(1).Color = ln2.Color;

% l = legend('NumColumns', 2);
l.Location = 'northwest';

 
xlabel('time to ripple center [ms]')

%%

plot_win_full = fs * 1.5; % +/- from events center, may want to zoom in more
plot_range_full = plot_center-plot_win_full:plot_center+plot_win_full;
plot_win_full = round(((plot_range_full - plot_center) / fs) * 1000);

subtightplot(7,4, 17+shift:4:24+shift, [0.04, 0.06])
yyaxis left
p = plot(plot_win_full, events.raw{ch}(rip_num, plot_range_full), 'k', 'LineWidth', 2.0);
hold on
ylim([min(events.raw{ch}(rip_num, plot_range_full)) max(events.raw{ch}(rip_num, plot_range_full))])

lowerleft = [plot_win(1), min(events.raw{ch}(rip_num, plot_range))];
width = 2*plot_win(end);

height = max(events.raw{ch}(rip_num, plot_range)) -  min(events.raw{ch}(rip_num, plot_range));
pos = [lowerleft(1), lowerleft(2), width, height];
r = rectangle('Position',pos);
r.EdgeColor = 0.5*[0 0.7 0];
r.LineWidth = 2.0;

% inds = find(events.goodRipples{ch} ~= events.goodRipples{ch}(rip_num));
% check = ((events.goodRipples{ch}(inds) > events.goodRipples{ch}(rip_num) - plot_win_full) & (events.goodRipples{ch}(inds) < events.goodRipples{ch}(rip_num) + plot_win_full));
% if any(check)
%     closeEvents = events.goodRipples{ch}(inds(check));
%     RelTime = closeEvents - events.goodRipples{ch}(rip_num);  
%     RelTimeRip = RelTime + plot_center;
%     p = plot(RelTime,events.raw{ch}(rip_num, RelTimeRip),'rs'); hold on;
%     p.MarkerFaceColor = 'r';
% end

ylabel('LFP \muV')    
vline(0)

yyaxis right

ln3 = plot(plot_win_full, abs(hilbert(events.rippleband{ch}(rip_num,plot_range_full))), 'g-'); hold on;
ln3.LineWidth = 1.5;
ln3.Color = [0 0.7 0];
ylabel('rippleband amp \muV')

% legend(ln3, {sprintf('%i-%iHz zscore', rippleband)})

% ylim([-5 5])
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = [0 .7 0];

% ax = gca;
% ax.XTickLabel = round(ax.XTick * (1/fs) * 1000); %xaxis in ms 
xlabel('time to ripple center [ms]')

fig = gcf;
fig.Color = [1 1 1];

return



