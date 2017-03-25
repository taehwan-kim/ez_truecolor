clc; close all;

dataforsave_450 = load('par_N_data_450.mat');
dataforsave_500 = load('par_N_data_500.mat');
dataforsave_550 = load('par_N_data_550.mat');
dataforsave_600 = load('par_N_data_600.mat');
dataforsave_450 = dataforsave_450.dataforsave;
dataforsave_500 = dataforsave_500.dataforsave;
dataforsave_550 = dataforsave_550.dataforsave;
dataforsave_600 = dataforsave_600.dataforsave;

dataformag = load('magnification_600.mat');
dataformag = dataformag.svfile;

dataforsave_z_450 = load('par_z_data_450.mat');
dataforsave_z_450 = dataforsave_z_450.dataforsave;
dataforsave_z_500 = load('par_z_data_500.mat');
dataforsave_z_500 = dataforsave_z_500.dataforsave;
dataforsave_z_550 = load('par_z_data_550.mat');
dataforsave_z_550 = dataforsave_z_550.dataforsave;
dataforsave_z_600 = load('par_z_data_600.mat');
dataforsave_z_600 = dataforsave_z_600.dataforsave;

dataforsave_bp_600 = load('par_z_bp_600.mat');
dataforsave_bp_600 = dataforsave_bp_600.dataforsave;

figure();

plot(dataforsave_450(1,:),1e9*dataforsave_450(2,:),'LineWidth',2);
hold on;
plot(dataforsave_500(1,:),1e9*dataforsave_500(2,:),'LineWidth',2);
hold on;
plot(dataforsave_550(1,:),1e9*dataforsave_550(2,:),'LineWidth',2);
hold on;
plot(dataforsave_600(1,:),1e9*dataforsave_600(2,:),'LineWidth',2);

set(gca, 'Xtick', 0:500:2000);
set(gca, 'Ytick', 0:2:10);
tk = get(gca,'XTick');
set(gca, 'Fontsize', 12);


title('CRLB_x (NA=1.4, Pixel size=100nm, Background=10)');
xlabel('Number of photons');
ylabel('Localization accuracy (nm)');

ax = gca;
set(ax.Title, 'Fontsize', 16);
set(ax.Title, 'FontWeight', 'Normal');
set(ax.XLabel, 'Fontsize', 14);
set(ax.YLabel, 'Fontsize', 14);

xlim([0 2200]);
ylim([0 10]);

lgd = legend('450nm','500nm','550nm','600nm');
lgd.Box = 'off';
lgd.FontSize = 14;


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, 'SaveFilex', '-dpdf');

figure();

plot(dataforsave_450(1,:),1e9*dataforsave_450(3,:),'LineWidth',2);
hold on;
plot(dataforsave_500(1,:),1e9*dataforsave_500(3,:),'LineWidth',2);
hold on;
plot(dataforsave_550(1,:),1e9*dataforsave_550(3,:),'LineWidth',2);
hold on;
plot(dataforsave_600(1,:),1e9*dataforsave_600(3,:),'LineWidth',2);

set(gca, 'Xtick', 0:500:2000);
set(gca, 'Ytick', 0:5:25);
tk = get(gca,'XTick');
set(gca, 'Fontsize', 12);

title('CRLB_\lambda (NA=1.4, Pixel size=100nm, Background=10)');
xlabel('Number of photons');
ylabel('Localization accuracy (nm)');

ax = gca;
set(ax.Title, 'Fontsize', 16);
set(ax.Title, 'FontWeight', 'Normal');
set(ax.XLabel, 'Fontsize', 14);
set(ax.YLabel, 'Fontsize', 14);
% 
% tk = get(gca,'XTick');
% set(gca, 'Fontsize', 12);
% tk = get(gca,'YTick');
% set(gca, 'Fontsize', 12);

xlim([0 2200]);
ylim([0 25]);

lgd = legend('450nm','500nm','550nm','600nm');
lgd.Box = 'off';
lgd.FontSize = 14;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, 'SaveFilel', '-dpdf');



figure();

h = plot(dataformag(1,:),dataformag(2,:),'LineWidth',2);

set(gca, 'Xtick', 0:50:500);
set(gca, 'Ytick', 1:0.2:2);
tk = get(gca,'XTick');
set(gca, 'Fontsize', 12);

title('PSF size dependency on axial position');
xlabel('\Deltaz (nm)');
ylabel('PSF Magnification Factor');

ax = gca;
set(ax.Title, 'Fontsize', 16);
set(ax.Title, 'FontWeight', 'Normal');
set(ax.XLabel, 'Fontsize', 14);
set(ax.YLabel, 'Fontsize', 14);

xlim([0 500]);
ylim([0.9 2]);

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, 'SaveFilemag', '-dpdf');



figure();

% plot(1e9*dataforsave_z_450(1,:),1e9*dataforsave_z_450(3,:),'LineWidth',2);
% hold on;
% plot(1e9*dataforsave_z_500(1,:),1e9*dataforsave_z_500(3,:),'LineWidth',2);
% hold on;
% plot(1e9*dataforsave_z_550(1,:),1e9*dataforsave_z_550(3,:),'LineWidth',2);
% hold on;
plot(1e9*dataforsave_z_600(1,:),1e9*dataforsave_z_600(3,:),'LineWidth',2);


set(gca, 'Xtick', 0:50:500);
set(gca, 'Ytick', 0:2:10);
tk = get(gca,'XTick');
set(gca, 'Fontsize', 12);

title('1/I(\lambda) (NA=1.4, Pixel size=100nm, Background=10)');
xlabel('\Deltaz (nm)');
ylabel('Localization accuracy (nm)');

ax = gca;
set(ax.Title, 'Fontsize', 16);
set(ax.Title, 'FontWeight', 'Normal');
set(ax.XLabel, 'Fontsize', 14);
set(ax.YLabel, 'Fontsize', 14);

xlim([0 500]);
ylim([0 10]);

% lgd = legend('450nm','500nm','550nm','600nm');
% lgd.Box = 'off';
% lgd.FontSize = 14;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, 'SaveFilez', '-dpdf');


figure();

% plot(1e9*dataforsave_z_450(1,:),1e9*dataforsave_z_450(3,:),'LineWidth',2);
% hold on;
% plot(1e9*dataforsave_z_500(1,:),1e9*dataforsave_z_500(3,:),'LineWidth',2);
% hold on;
% plot(1e9*dataforsave_z_550(1,:),1e9*dataforsave_z_550(3,:),'LineWidth',2);
% hold on;
plot(1e9*[-1*fliplr(dataforsave_bp_600(1,2:end)) dataforsave_bp_600(1,:)],1e9*[fliplr(dataforsave_bp_600(3,2:end)) dataforsave_bp_600(3,:)],'LineWidth',2);


set(gca, 'Xtick', -500:100:500);
set(gca, 'Ytick', 0:5:15);
tk = get(gca,'XTick');
set(gca, 'Fontsize', 12);

title('CRLB_\lambda (NA=1.4, Pixel size=100nm, Background=10)');
xlabel('\Deltaz (nm)');
ylabel('Localization accuracy (nm)');

ax = gca;
set(ax.Title, 'Fontsize', 16);
set(ax.Title, 'FontWeight', 'Normal');
set(ax.XLabel, 'Fontsize', 14);
set(ax.YLabel, 'Fontsize', 14);

xlim([-500 500]);
ylim([0 15]);

% lgd = legend('450nm','500nm','550nm','600nm');
% lgd.Box = 'off';
% lgd.FontSize = 14;

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig, 'SaveFilebp', '-dpdf');


% zlocations = dataforsave_z(1,:);
% crlb = dataforsave_z(3,:).^2;
% lambda = 600e-9;
% haha = dataformag(2,1:11).^2.*crlb + lambda^2 * (dataformag(2,1:11)-1).^2;
