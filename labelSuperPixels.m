clear all
close all
clc

I1=imread('sec.jpg');
I2=imread('label.jpg');
[h,w,s]=size(I1);
   for i=1:h
      for j=1:w
          if I2(i,j,1)==0
              I1(i,j,1)=0;
              I1(i,j,2)=0;
              I1(i,j,3)=255;
          end     
      end
   end

figure
imshow(I1)