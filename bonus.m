% Exercice Bonus
% workspace preparation
clear all
clc
close all

% Bonus
% Récupération des données
% Bonus
% Récupération des données
% a1=[0,00000000e+00
% 1,66656608e-01
% 3,00000000e-01
% 6,00000000e-01
% 9,00000000e-01
% 1,20000000e+00
% 1,50000000e+00
% ]
% b1=[2,98000000e+02
% 3,98000000e+02
% 5,34749046e+02
% 1,90090090e+03
% 1,96599825e+03
% 1,97297297e+03
% 1,97297297e+03
% ]
% plot(a1,b1)

img1 = imread('1.png');
img2 = imread('2.png');
img3 = imread('3.png');
img4 = imread('4.png');

myVideo = VideoWriter('myvideo.avi');
myVideo.FrameRate = 30;
open(myVideo);
% % resize the image to 914 by 744
% img_resized = imresize(img2, [744 914]);

% % save the resized image
% imwrite(img_resized, 'img2re.jpg');
writeVideo(myVideo, img1);
writeVideo(myVideo, img2);
writeVideo(myVideo, img3);
writeVideo(myVideo, img4);

close(myVideo);
