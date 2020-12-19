function klabels = PerformSuperpixelSLIC(img,img_Lab, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, STEP, compactness)

[m_height, m_width, m_channel] = size(img_Lab);
numseeds = size(kseedsl);
img_Lab = double(img_Lab);
klabels = zeros(m_height, m_width);
degree=3;
%-----------------Amap, a feature used for SLIC----------------------------
Amap=zeros(m_height, m_width);
kh=floor(m_height/degree);
kw=floor(m_width/degree);

 new1=img(:,:,1);
 new2=img(:,:,2);
 new3=img(:,:,3);
 for op1=1:degree:(kh-1)*degree+1
    for op2=1:degree:(kw-1)*degree+1
        I_1(:,:,1)=double(new1(op1:op1+degree-1,op2:op2+degree-1));
        I_1(:,:,2)=double(new2(op1:op1+degree-1,op2:op2+degree-1));
        I_1(:,:,3)=double(new3(op1:op1+degree-1,op2:op2+degree-1));
        I_1=uint8(I_1);
         for i=1:degree
            for j=1:degree
                dark_I(i,j)=min(I_1(i,j,:));
            end;
        end;
        dark_col=reshape(dark_I,degree*degree,1);
        I_col=reshape(I_1,degree*degree,3);
        [dark_sort,Index]=sort(dark_col,'descend');
        num_A=floor(degree*degree/2);%%%%%patch_size>10，divided by 100；4<patch_size<10，divided /10
        Index_A=Index(1:num_A);
        A=max(max(I_col(Index_A,:)));
        A=double(A);
        for ax=op1:op1+degree-1
            for ay=op2:op2+degree-1
                Amap(ax,ay)=A;
            end
        end
    end
 end
%--------------------------------------------------------

invwt = 1/((double(STEP)/double(compactness))*(double(STEP)/double(compactness)));
numk = numseeds;
for itr = 1: 10   
    sigmal = zeros(numseeds(1,1), 1);
    sigmaa = zeros(numseeds(1,1), 1);
    sigmab = zeros(numseeds(1,1), 1);
    sigmax = zeros(numseeds(1,1), 1);
    sigmay = zeros(numseeds(1,1), 1);
    clustersize = zeros(numseeds(1,1), 1);
    inv = zeros(numseeds(1,1), 1);
    distvec = double(100000*ones(m_height, m_width));
   
    for n = 1: numk
        y1 = max(1, kseedsy(n, 1)-STEP);
        y2 = min(m_height, kseedsy(n, 1)+STEP);
        x1 = max(1, kseedsx(n, 1)-STEP);
        x2 = min(m_width, kseedsx(n, 1)+STEP);
       
        for y = y1: y2
            for x = x1: x2
                %dist_lab = abs(img_Lab(y, x, 1)-kseedsl(n))+abs(img_Lab(y, x, 2)-kseedsa(n))+abs(img_Lab(y, x, 3)-kseedsb(n));
                dist_lab = (img_Lab(y, x, 1)-kseedsl(n, 1))^2+(img_Lab(y, x, 2)-kseedsa(n, 1))^2+(img_Lab(y, x, 3)-kseedsb(n, 1))^2;
                dist_xy = (double(y)-kseedsy(n, 1))*(double(y)-kseedsy(n, 1)) + (double(x)-kseedsx(n, 1))*(double(x)-kseedsx(n, 1));
%                  Da=abs(Amap(y,x)-Amap(fix(kseedsy(n, 1)),fix(kseedsx(n, 1))));
                Da=(Amap(y,x)-Amap(fix(kseedsy(n, 1)),fix(kseedsx(n, 1))))^2;
                dist = dist_lab + dist_xy*invwt+Da;
                
                %m = (y-1)*m_width+x;
                if (dist<distvec(y, x))
                    distvec(y, x) = dist;
                    klabels(y, x) = n;
                end
            end
        end
    end
    
    ind = 1;
    for r = 1: m_height
        for c = 1: m_width
            sigmal(klabels(r, c),1) = sigmal(klabels(r, c),1)+img_Lab(r, c, 1);
            sigmaa(klabels(r, c),1) = sigmaa(klabels(r, c),1)+img_Lab(r, c, 2);
            sigmab(klabels(r, c),1) = sigmab(klabels(r, c),1)+img_Lab(r, c, 3);
            sigmax(klabels(r, c),1) = sigmax(klabels(r, c),1)+c;
            sigmay(klabels(r, c),1) = sigmay(klabels(r, c),1)+r;
            clustersize(klabels(r, c),1) = clustersize(klabels(r, c),1)+1;
        end
    end
    for m = 1: numseeds
        if (clustersize(m, 1)<=0)
            clustersize(m, 1) = 1;
        end
        inv(m, 1) = 1/clustersize(m, 1);
    end
    for m = 1: numseeds
        kseedsl(m, 1) = sigmal(m, 1)*inv(m, 1);
        kseedsa(m, 1) = sigmaa(m, 1)*inv(m, 1);
        kseedsb(m, 1) = sigmab(m, 1)*inv(m, 1);
        kseedsx(m, 1) = sigmax(m, 1)*inv(m, 1);
        kseedsy(m, 1) = sigmay(m, 1)*inv(m, 1);
    end
end
end


