clear all; 
close all;

temp = load('width_03_10.mat');
width_03_10 = temp.width;
temp = load('width_03_15.mat');
width_03_15 = temp.width;
temp = load('inten_03_10.mat');
inten_03_10 = temp.inten;
temp = load('inten_03_15.mat');
inten_03_15 = temp.inten;

figure(1);
subplot(2,1,1);
hist(width_03_10);
xlim([0.9 0.95]);
subplot(2,1,2);
histfit(width_03_15);
xlim([0.9 0.95]);

figure(2);
subplot(2,1,1);
histfit(inten_03_10);
xlim([0.25 0.35]);
subplot(2,1,2);
histfit(inten_03_15);
xlim([0.25 0.35]);
