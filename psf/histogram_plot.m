clear all; 
close all;

temp = load('width_005_10.mat');
width_005_10 = temp.width;
temp = load('inten_005_10.mat');
inten_005_10 = temp.inten;
temp = load('width_005_20.mat');
width_005_20 = temp.width;
temp = load('inten_010_10.mat');
inten_010_10 = temp.inten;

figure(1);
subplot(2,1,1);
histfit(width_005_10);
% xlim([0.9 1.3]);
subplot(2,1,2);
histfit(width_005_20);
% xlim([0.9 1.3]);

figure(2);
subplot(2,1,1);
histfit(inten_005_10);
% xlim([0.25 0.35]);
subplot(2,1,2);
histfit(inten_010_10);
% xlim([0.25 0.35]);
