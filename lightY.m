clc;
clear all
close all

LR=imread('outdoor/45-11.jpg');
GT=imread('outdoor/45-11.jpg');
%figure
%imshow(LR);
im_l_ycbcr = rgb2ycbcr(LR);
im_l_y = im_l_ycbcr(:, :, 1);
im_l_cb = im_l_ycbcr(:, :, 2);
im_l_cr = im_l_ycbcr(:, :, 3);
im_l_y=double(im_l_y);
im_l_cb=double(im_l_cb);
im_l_cr=double(im_l_cr);
[h,w]=size(im_l_y);
m=0;
for i=1:h
    for j=1:w
        m=m+im_l_y(i,j);
    end
end
n=m/(h*w);
disp(n);
im_2_ycbcr = rgb2ycbcr(GT);
im_2_y = im_2_ycbcr(:, :, 1);
im_2_cb = im_2_ycbcr(:, :, 2);
im_2_cr = im_2_ycbcr(:, :, 3);
im_2_y=double(im_2_y);
im_2_cb=double(im_2_cb);
im_2_cr=double(im_2_cr);
m1=0;
for i=1:h
    for j=1:w
        m1=m1+im_2_y(i,j);
    end
end
n1=m1/(h*w);
disp(n1);
% p=n1/n;
p=1.8;
disp(p)
im_h_y=im_l_y*p-80;
im_h_ycbcr(:, :, 1) = im_h_y;
im_h_ycbcr(:, :, 2) = im_l_cb;
im_h_ycbcr(:, :, 3) = im_l_cr;
im_h = ycbcr2rgb(uint8(im_h_ycbcr));
figure
imshow(im_h);
imwrite(im_h,'35-111.jpg');

