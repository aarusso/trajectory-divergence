%% DIVERGENCE DEMO
%
% Data formatting
% For tangling and divergence analyses, data should be formatted as for
% jPCA. Each data structure (e.g. D_m1, loacted in M1_sampleData.mat)
% contains C elements corresponding to the number fo conditions (here, 2).
% The t x n matrix 'A' contains the trial-averaged firing rates as a
% function of time for each neuron. All other fields are optional. The
% 'times' field indicates the time course of each sample in A and the
% 'analyzeTimes' field indicates which of these times should be analyzed
% (e.g. by jPCA, tangling, or divergence).
% 
% Sample data contents The data included in this demo were recorded in the
% Churchland lab at Columbia University. Single-unit neural data were
% recorded sequentially from a macaque monkey during the performance of a
% hand-pedaling task. Monkeys grasped a hand crank and cycled through a
% virtual environment for a number of prescribed cycles as indicated by a
% visual cue. The data included here correspond to two seven-cycle
% conditions (forward bottom-start and backward bottom-start for monkey C).
%  For more information about data processing and the task, see Russo et al
%  Neuron 2020.
%
%
% AUTHORSHIP
% |Author: Abigail Russo,|
% |Email: russoaa@gmail.com,|
% |Dated: May 2019|
%% Plot example behavioral kinamtics
% The sample data contain two conditions: seven-cycle for forward pedaling
% starting at the bottom and seven cycle for backward pedaling starting at
% the bottom. Plot the kinematic data (hand position and velocity, torque,
% and wrist angle) from 8 example sessions to get a feel for the behavior.
% Data have been trial averaged and aligned as described in Russo et al
% Neuron 2020. The analysis time window is 100ms before movement onset to
% 100ms after movement offset.

load kinematic_sampleData


T = length(p(1).A);
ticks = sort([500 1500:500:T T-500]);
tickLabels = arrayfun(@(n) {num2str(n)}, (ticks/1000) - 1.5); % time from movement onset in seconds
tickLabels{1} = 'target on';
tickLabels{end-1} = 'stop';

figure; 
subplot(421); hold on;
h = plot(p(1).A,'color',[.2 .8 .2]);
j = plot(p(2).A,'color',[.8 .2 .2]);
ylabel('world position')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])
analyzeTimes = p(1).analyzeTimes;
k = plot([analyzeTimes(1) analyzeTimes(end)],[0 0],'k');
plot([analyzeTimes(1) analyzeTimes(1)],[-1 1],'k');
plot([analyzeTimes(end) analyzeTimes(end)],[-1 1],'k');
legend([h(1) j(1) k],{'Forward bottom start','Backward bottom start','Analysis time window'},'Location','NW')

subplot(423); hold on;
plot(diff(p(1).A),'color',[.2 .8 .2])
plot(diff(p(2).A),'color',[.8 .2 .2])
ylabel('angular hand velocity')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(425); hold on;
plot(tq(1).A,'color',[.2 .8 .2])
plot(tq(2).A,'color',[.8 .2 .2])
ylabel('torque')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(427); hold on;
plot(w(1).A,'color',[.2 .8 .2])
plot(w(2).A,'color',[.8 .2 .2])
ylabel('wrist angle')
set(gca,'XTick',ticks,'XTickLabels',tickLabels,'XLim',[0 T])
xlabel('time from movement onset (s)')


subplot(422); hold on;
plot(vp(1).A,'color',[.2 .8 .2])
plot(vp(2).A,'color',[.8 .2 .2])
ylabel('vertical hand position')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(424); hold on;
plot(hp(1).A,'color',[.2 .8 .2])
plot(hp(2).A,'color',[.8 .2 .2])
ylabel('horizontal hand position')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(426); hold on;
plot(vv(1).A,'color',[.2 .8 .2])
plot(vv(2).A,'color',[.8 .2 .2])
ylabel('vertical hand velocity')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(428); hold on;
plot(hv(1).A,'color',[.2 .8 .2])
plot(hv(2).A,'color',[.8 .2 .2])
ylabel('horizontal hand velocity')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])
set(gca,'XTick',ticks,'XTickLabels',tickLabels)
xlabel('time from movement onset (s)')

sgtitle('Example kinematics from 8 recording sessions')


%% Analyze divergence
load M1_sampleData
load SMA_sampleData



Dvg_m1 = divergeAnalysis( D_m1, 'softenNorm',5);
Dvg_sma = divergeAnalysis(D_sma,'softenNorm',5);

% make scatter plot with sample data, formatted as in Russo et al 2020 Fig 7A
formattedScatter(Dvg_sma, Dvg_m1, {'SMA Divergence','M1 Divergence'});


% make histogram difference plot with sample data, formatted as in Russo et al 2020 Fig 7C
figure; hold on
Dvg_Diff = Dvg_sma-Dvg_m1;
myRange = [min(Dvg_Diff) max(Dvg_Diff)];
myEdge = max(abs(myRange)); bw = 10;
myEdge = ceil(myEdge/10)*10; 
histParams = {'DisplayStyle','bar','FaceColor',[.2 .5 .8],...
   'EdgeAlpha',0,'BinEdges',[-myEdge:bw:myEdge],'Normalization','probability'};
histogram(Dvg_Diff,histParams{:})
xlim = max(abs(get(gca,'XLim')));
set(gca,'XLim',[-xlim xlim])
ylim = get(gca,'YLim');
plot([0 0],[0 ylim(2)],'--','Color','k','LineWidth',2)

params = struct();
params.fontSize = 14;
params.axisLabel = 'Divergence (SMA - M1)';
params.axisOffset = 0;
AxisMMC(-myEdge, myEdge, params);

params.axisLabel = 'Fraction of time points';
params.axisOrientation = 'v';
params.axisOffset = -xlim;

AxisMMC(0, ylim(2), params);
set(gca,'Visible','off');
set(gcf,'Color','w')
axis([-xlim-bw, xlim+bw, -.01, ylim(2)*1.05])