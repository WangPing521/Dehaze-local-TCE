function [light,Anum,img_ContoursEX] = DrawContoursAroundSegments_EX(img, klabels,K)

[h,w]=size(klabels);
Anum=zeros(K+1,1);
light=zeros(K+1,1);
degree=3;
ww=0.25;
t0=0.1;
m=0.3;
%------------------------------------------------
dx1 = [-1, 0, 1, -1, 1, -1, 0, 1];
dy2 = [-1, -1, -1, 0, 0, 1, 1, 1];
img_ContoursEX1 = img;
for j = 1: h
    for k = 1: w
        np = 0;
        for i = 1: 8
            x = k+dx1(1, i);
            y = j+dy2(1, i);
            if (x>0&&x<=w&&y>0&&y<h)
               if (klabels(j, k)~=klabels(y, x))
                   np = np+1;
               end
            end
        end
        if (np>2)
            img_ContoursEX1(j, k, 1) = 255;
            img_ContoursEX1(j, k, 2) = 255;
            img_ContoursEX1(j, k, 3) = 255;
        end
    end
end
figure
imshow(img_ContoursEX1)
title('SLIC results');
%求解每个超像素区域内的局部大气光值
for sg=1:K-1
    n=1;
    for pi=1:h
        for pj=1:w
            if(klabels(pi,pj)==sg)
                newimg(n,1,:)=img(pi,pj,:);
                n=n+1;
            end
        end
    end
    for ai=1:n-1
      dark_img(ai,1)=min(newimg(ai,1,:));
    end
    [dark_sort,Index]=sort(dark_img,'descend');
    num_A=floor((n-1)/10);
    Index_A=Index(1:num_A);
    A=max(max(newimg(Index_A,:)));
    A=double(A);
    Anum(sg,1)=A;
%求解每个超像素的平均亮度（作为超像素合并的准则）
    new=newimg;
    newgray=rgb2gray(new);
    newgray=double(newgray);
    sum=0;
    for sum1=1:n-1
        sum=sum+newgray(sum1,1);
    end
    light(sg,1)=sum/(n-1);
end

lightmax=max(light);
lightmin=min(light);
distance=fix((lightmax-lightmin)/5);
% fprintf('%3d',distance);
%------------------超像素合并，并计算合并后每个区域内的大气光值-------------------------
dx = [-1, 0, 1, -1, 1, -1, 0, 1];
dy = [-1, -1, -1, 0, 0, 1, 1, 1];
    for mr1=2:h-1
        for mr2=2:w-1
            for mri = 1: 8
                y1 = mr1+dy(1, mri);
                x1 = mr2+dx(1, mri);
                if (x1>0&&x1<=w&&y1>0&&y1<h)
                   if (klabels(mr1, mr2)~=klabels(y1, x1))
                       a=klabels(mr1, mr2);
                       b=klabels(y1, x1);
                       if(abs(light(a,1)-light(b,1))<distance)
                           %遍历，将标签为b的像素的标签都置为a
                           for cor1=1:h
                               for cor2=1:w
                                   if(klabels(cor1,cor2)==b)
                                       klabels(cor1,cor2)=a;
                                   end
                               end
                           end                       
                           Anum(a,1)=max(Anum(a,1),Anum(b,1));
                           Anum(b,1)=Anum(a,1);
                       end
                   end
                end
            end        
        end
    end
% ==========去雾函数（写进ThreeAirLightHaze里）=================
tic
out=ThreeAirLightHaze1(img,degree,m,ww,t0,klabels,Anum);
toc
%-------------------------------------------------
img_ContoursEX = img;
for j = 1: h
    for k = 1: w
        np = 0;
        for i = 1: 8
            x = k+dx(1, i);
            y = j+dy(1, i);
            if (x>0&&x<=w&&y>0&&y<h)
               if (klabels(j, k)~=klabels(y, x))
                   np = np+1;
               end
            end
        end
        if (np>2)
            img_ContoursEX(j, k, 1) = 0;
            img_ContoursEX(j, k, 2) = 0;
            img_ContoursEX(j, k, 3) = 255;
        end
    end
end

% I1=imread('sec.jpg');
% [h,w,s]=size(I1);
%    for i=1:h
%       for j=1:w
%           if img_ContoursEX(i,j,1)==0
%               I1(i,j,1)=0;
%               I1(i,j,2)=0;
%               I1(i,j,3)=255;
%           end     
%       end
%    end


figure;
imshow(out)
title('dehazed image')
imwrite(out,'outdoor\45-11.jpg');

% figure;
% imshow(img_ContoursEX)
% imwrite(img_ContoursEX,'label.jpg');
% title('superpixels')
%------------------------------------------------------
