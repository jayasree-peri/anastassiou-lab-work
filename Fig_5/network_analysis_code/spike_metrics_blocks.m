
%% load spike groups to be analysed

% WG1
composition_v = {'10001' '10002' '10003' '10004' '10005'};

% WG4
%composition_v = {'1001' '1002' '1003' '1004' '1005'};

% WG4 x2syn
%composition_v = {'1021' '1022' '1023' '1024' '1025'};

% WG4 alt
%composition_v = {'1011' '1012' '1013' '1014' '1015'};

%%


%% set up parameters and analyze excitatory cells

% define a few paths and some parameters
path_s = cell(length(composition_v),1);
for kk = 1:length(path_s)
    path_s{kk} = ['./' num2str(composition_v{kk}) '/'];
    compo_v(kk) = str2num(composition_v{kk});
end
spike_f = 'spikes.csv';

% some more parameters...
bin_dt = 10;                    % binning time step for assessing spike frequency
pre = 0;
sigma_times = 1.5;              % sigma for defining the event detection threshold
dt_event_detection = 200;       % time interval to start detection of next event
flag_plot = 0;                  % flag to plot event detection
flag_plot2 = 0;                 % flag to plot event maxima

max_detect_win =  1;             % number of points after threshold crossing to look for max
max_win = 12;                    % number of bins around the max to calculate stats

% number of cell types to be considered
CT = cell(2,1);
CT{1} = [1 499];
CT{2} = [500 505];

% Extract some data for analysis...
thresh = zeros(length(path_s), 1);
bins = cell(length(composition_v),1);
f_sp = cell(length(composition_v),1);
f_sp_ct = cell(length(composition_v),1);
ind_event_start = cell(length(composition_v),1);
ind_event_max = cell(length(composition_v),1);
win_max_left = cell(length(composition_v),1);
win_max_right = cell(length(composition_v),1);
fit_param = cell(length(composition_v),1);

% analyze the propeprties of excitatory cells


for kk = 1:length(path_s)
    [post, spikeTimes, spiked_v, ct_ind, bins{kk}, f_sp_ct{kk}, f_sp{kk}] = spike_metrics_industrial(path_s{kk}, spike_f, bin_dt, pre, 1, CT);  
    thresh(kk) = mean(f_sp{kk}) + sigma_times.*std(f_sp{kk});
    % detect threshold crossing...
    [ind] = detect_threshold_crossing(kk, thresh, dt_event_detection, bin_dt, bins, f_sp, flag_plot);
    ind_event_start{kk}(1:length(ind)) = ind;
    % detect maximum fsp at every event
    [ind_m] = find_max(kk, ind_event_start, max_detect_win, bins, f_sp, flag_plot2);
    ind_event_max{kk}(1:length(ind_m)) = ind_m;
    % determine window (bins) around maximum to fit function (gaussian)
    [win_left, win_right] = bin_win_fit(kk, max_win, bins, ind_event_max{kk});
    win_max_left{kk}(1:length(win_left)) = win_left;
    win_max_right{kk}(1:length(win_right)) = win_right;
    % fit function at every event and extract some numbers...
    [fit_params] = fit_data_func(kk, win_max_left, win_max_right, bins, f_sp, 'gauss1');
    fit_param{kk}(1:size(fit_params,1), 1:size(fit_params,2)) = fit_params;
end



%% analyze the propeprties of inhibitory cells

% rename the variables
for i = 1:length(f_sp)
    f_sp{i}=f_sp_ct{i}(:,2);
end

for kk = 1:length(path_s)
%    [post, spikeTimes, spiked_v, ct_ind, bins{kk}, f_sp_ct{kk}, f_sp{kk}] = spike_metrics_industrial(path_s{kk}, spike_f, bin_dt, pre, 1, CT);  
    thresh(kk) = mean(f_sp{kk}) + sigma_times.*std(f_sp{kk});
    % detect threshold crossing...
    [ind] = detect_threshold_crossing(kk, thresh, dt_event_detection, bin_dt, bins, f_sp, flag_plot);
    ind_event_start{kk}(1:length(ind)) = ind;
    % detect maximum fsp at every event
    [ind_m] = find_max(kk, ind_event_start, max_detect_win, bins, f_sp, flag_plot2);
    ind_event_max{kk}(1:length(ind_m)) = ind_m;
    % determine window (bins) around maximum to fit function (gaussian)
    [win_left, win_right] = bin_win_fit(kk, max_win, bins, ind_event_max{kk});
    win_max_left{kk}(1:length(win_left)) = win_left;
    win_max_right{kk}(1:length(win_right)) = win_right;
    % fit function at every event and extract some numbers...
    [fit_params] = fit_data_func(kk, win_max_left, win_max_right, bins, f_sp, 'gauss1');
    fit_param{kk}(1:size(fit_params,1), 1:size(fit_params,2)) = fit_params;
end
%%



%% analyze the relative timing between the exc and inh populations

delay_exc_inh=length(fit_params_exc);

delay_threshold=5;

for i=1:1:length(fit_params_exc)
   
   delay=fit_params_exc(i,2) - fit_params_inh(i,2);

   % compute the delay if there is no difference
   if abs(delay) <= delay_threshold
   delay_exc_inh(i)=fit_params_exc(i,2) - fit_params_inh(i,2);   
   end
               
   % compute the delay if there is too large difference (different events for exc and inh cells)
   if abs(delay) >= delay_threshold
       
       delay_1=fit_params_exc(i,2) - fit_params_inh(i+1,2);       
       delay_exc_inh(i)=fit_params_exc(i,2) - fit_params_inh(i+1,2);
       
       if abs(delay_1) >= delay_threshold
        delay_2=fit_params_exc(i,2) - fit_params_inh(i+2,2);       
        delay_exc_inh(i)=fit_params_exc(i,2) - fit_params_inh(i+2,2);
       end
       
   end
           
end


boxplot(delay_exc_inh);
set(gca,'Fontsize',20);

box off;
ylabel('Exc-to-Inh delay (ms)');
ylim([-1,2]);
xticklabels('WG4 2syn')


%%




%% plot events statistics



%% Plot the raster and set up figure parameters

linew = 2.0;
ax_linew = 1.5;
msize = 12.0;
ax_msize = 15.0;
color_v = ['b' 'r'];

figure(1000); set(gcf,'color','w');

for kk = 1:length(path_s)
    
    [post, spikeTimes, spiked_v, ct_ind] = spike_metrics_industrial(path_s{kk}, spike_f, bin_dt, pre, 1, CT);
        
    % raster plot...
    subplot(length(path_s), 2, 2*(kk-1)+1);
    for j=1:size(ct_ind,1)
        for i = 1:1:length(ct_ind{j})
            plot(round(spikeTimes{ct_ind{j}(i)}) + pre, spiked_v(ct_ind{j}(i)),'.', 'color', color_v(j), 'MarkerSize',8);
            %xlim([1000 5000]);
            hold on;
        end
    end
    
    xlabel('time / ms','fontsize', 25);
    ylabel('Cell ID','fontsize', 25);
    box off; set(gca,'linewidth',ax_linew);
    axis([-pre-25 post+25 -25 spiked_v(length(spiked_v))+25]);
    
    subplot(length(path_s), 2, 2*(kk-1)+2);
    for j=1:size(ct_ind,1)
        plot(bins{kk}, f_sp_ct{kk}(:,j), '-', 'color', color_v(j), 'linewidth', linew);
        %xlim([1000 5000]);
        hold on;
    end
    plot(bins{kk}, f_sp{kk}, '-k', 'linewidth', linew-1.5);
    plot([-pre-25 post+25], [thresh(kk) thresh(kk)], '--', 'color', [0.3 0.9 0.3], 'linewidth', linew-0.5);
    ylabel('Spike frequency / Hz','fontsize', 14);
    xlabel('time / ms','fontsize', 14);
    box off; set(gca,'linewidth',ax_linew);
    xlim([-pre-25 post+25]);
end

%%


%% Plot event statistics 

fit_a = zeros(length(fit_param),1); fit_o = zeros(length(fit_param),1); fit_w = zeros(length(fit_param),1);
for j = 1:length(fit_param), fit_a(j) = fit_param{j}(1); fit_o(j) = fit_param{j}(2);  fit_w(j) = fit_param{j}(3); end

figure(3043); set(gcf,'color','w');
subplot(1,3,1);
plot(compo_v, fit_a, '-ok', 'linewidth', linew); hold on;
box off;
%axis([-0.05 1.05 -5 105]);
ylabel('event spike frequency amplitude (Hz)','fontsize', 14);
xlabel('WG1 vs. WG4 composition','fontsize', 14);

subplot(1,3,2);
plot(compo_v, fit_o, '-ok', 'linewidth', linew); hold on;
box off;
ylabel('event centroid (timing) / ms','fontsize', 14);
xlabel('WG1 vs. WG4 composition','fontsize', 14);

subplot(1,3,3);
plot(compo_v, fit_w, '-ok', 'linewidth', linew); hold on;
box off;
ylabel('event duration / ms','fontsize', 14);
xlabel('WG1 vs. WG4 composition','fontsize', 14);

%%


%% compute the the mean half-width, amplitude and number of events

% assign fit param to the particular one

fit_param=fit_param;

for j=1:length(fit_param)
    
    % compute the means
    n_events(j)=length(fit_param{j}(:,1));
    mean_a(j)=mean(fit_param{j}(:,1));
    mean_fwth(j)=mean(fit_param{j}(:,3))*2*sqrt(2*log(2))*bin_dt;
    
end

% pooling the data

% mean and variance of the amplitude

a_100_2syn=mean_a;        % Hz, peak frequency
hz_100_2syn=n_events/10;  % Hz, event frequency
fwth_100_2syn=mean_fwth;  % ms, event duration

%a_mean_100=mean(mean_a);
%a_mean_std_100=std(mean_a);

% mean number of events 
%events_mean_1=mean(n_events);
%events_std_1=std(n_events);

% mean and std of width
%width_mean_1=mean(mean_fwth);
%width_std_1=std(mean_fwth);


%%


%% WG4 proportion box plots

figure
boxplot([a_0',a_8',a_17',a_25',a_33',a_42',a_50',a_58',a_67',a_75',a_83',a_92',a_100'],{'0','8','17','25','33','42','50','58','67','75','83','92','100'})
set(gca,'Fontsize',25)
xlabel('WG4 cells / %')
ylabel('Peak frequency / Hz')
box off

figure
boxplot([hz_0',hz_8',hz_17',hz_25',hz_33',hz_42',hz_50',hz_58',hz_67',hz_75',hz_83',hz_92',hz_100'],{'0','8','17','25','33','42','50','58','67','75','83','92','100'})
set(gca,'Fontsize',25)
xlabel('WG4 cells / %')
ylabel('Event frequency / Hz')
box off


figure
boxplot([fwth_0',fwth_8',fwth_17',fwth_25',fwth_33',fwth_42',fwth_50',fwth_58',fwth_67',fwth_75',fwth_83',fwth_92',fwth_100'],{'0','8','17','25','33','42','50','58','67','75','83','92','100'})
set(gca,'Fontsize',25)
xlabel('WG4 cells / %')
ylabel('Event duration / ms')
box off


%%


%% WG1 vs WG4 vs WG4_double reccurent spines for excitatory cells

figure
%subplot(3,1,1)
boxplot([a_0_2syn',a_100_2syn'],{'0','100'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG1', 'WG4'})
ylabel('Peak amplitude / Hz')
box off

figure
%subplot(3,1,2)
boxplot([hz_0_2syn',hz_100_2syn'],{'0','100'})
set(gca,'Fontsize',25)
xlabel('WG4 cells / %')
xticklabels({'WG1', 'WG4'})
ylabel('Event frequency / Hz')
box off

figure
%subplot(3,1,3)
boxplot([fwth_0_2syn',fwth_100_2syn'],{'0','100'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG1', 'WG4'})
ylabel('Event duration / ms')
box off


%%


%% Plot WG4 2x syn altered

figure
%subplot(3,1,1)
boxplot([hz_10081_alt'],{'100+altered'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG4 2syn altered'})
ylabel('Event frequency / Hz')
%ylim([1.5;4.5]);
box off


figure
%subplot(3,1,1)
boxplot([fwth_10081_alt'],{'100+altered'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG4 2syn altered'})
ylabel('Event duration / Hz')
%ylim([39.20;72]);
box off


figure
%subplot(3,1,1)
boxplot([a_10081_alt'],{'100+altered'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG4 2syn altered'})
ylabel('Event amplitude / Hz')
%ylim([20;80]);
box off


%%



%% WG1 vs WG4 vs WG4_double reccurent spines for inhibitory cells

figure
%subplot(3,1,1)
boxplot([a_0i',a_100i',a_1002i'],{'0','100','100+2rec'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG1', 'WG4' 'WG4 2syn'})
ylabel('Peak frequency / Hz')
box off

figure
%subplot(3,1,2)
boxplot([hz_0i',hz_100i',hz_1002i'],{'0','100','100+2rec'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG1', 'WG4' 'WG4 2syn'})
ylabel('Event frequency / Hz')
box off

figure
%subplot(3,1,3)
boxplot([fwth_0i',fwth_100i',fwth_1002i'],{'0','100','100+2rec'})
set(gca,'Fontsize',25)
%xlabel('WG4 cells / %')
xticklabels({'WG1', 'WG4' 'WG4 2syn'})
ylabel('Event duration / ms')
box off


%%