% dehaze
function [out]=ThreeAirLightHaze1(img,degree,m,ww,t0,klabels,Anum)

[h,w,s]=size(img);

for itr=1:5

   
    under_50=0;
    for i=1:h
        for j=1:w
            dark_I(i,j)=min(img(i,j,:));
        end;
    end;
    for i=1:h
        for j=1:w
            if(dark_I(i,j)<50)
                under_50=under_50+1;
            end;
        end;
    end;
    total=size(img,1)*size(img,2)*size(img,3); 
    percent=under_50/total;
    if(percent<m)          
        minV=255;
        maxV=0;
        countlabel=0;
        img=double(img);

        
        dark_I1=dark_I;
        kh=floor(h/degree);
        kw=floor(w/degree);
        lh=mod(h,degree);
        lw=mod(w,degree);

        for op1=1:degree:(kh-1)*degree+1
            for op2=1:degree:(kw-1)*degree+1
                darkI1_1=double(dark_I1(op1:op1+degree-1,op2:op2+degree-1));
                label=min(min(darkI1_1));
                for oop1=1:degree
                    for oop2=1:degree
                        darkI1_1(oop1,oop2)=label;
                    end;
                end;
                countlabel=countlabel+label;
                if label<minV
                    minV=label;
                end;
                if label>maxV
                    maxV=label;
                end;
                %splice
                if(op2==1)
                    I_linkcol=darkI1_1;
                else
                    I_linkcol=cat(2,I_linkcol,darkI1_1);
                end;
            end;
            if(op1==1)
                I_linkrow=I_linkcol;
            else
                I_linkrow=cat(1,I_linkrow,I_linkcol);
            end;
        end;
        I_row=I_linkrow;
        % splice
        if(lh==0&&lw==0)
            I_row=uint8(I_row);
        end;
        if(lh==0&&lw~=0)
            rightlink=double(dark_I(1:h,w-lw+1:w));
            I_row=cat(2,I_row,rightlink);
            I_row=uint8(I_row);
        end;
        if(lh~=0&&lw==0)
            downlink=double(dark_I(h-lh+1:h,1:w));
             I_row=cat(1,I_row,downlink);
             I_row=uint8(I_row);
        end;    
        if(lh~=0&&lw~=0)
            rightlink=double(dark_I(1:h-lh,w-lw+1:w));
            downlink=double(dark_I(h-lh+1:h,1:w));
            I_row=cat(2,I_row,rightlink);
            I_row=cat(1,I_row,downlink);
            I_row=uint8(I_row);
        end;
        % =======================================
         dark_channel=double(I_row);
         for opp1=1:degree:(kh-1)*degree+1
            for opp2=1:degree:(kw-1)*degree+1
                dark_channel1=double(dark_channel(opp1:opp1+degree-1,opp2:opp2+degree-1));
%                 numa=klabels(opp1,opp2);
%                 A=Anum(numa,1);
                  
                  A=150; 
                  t1=1-ww*(dark_channel1/A);
                %splice
                if(opp2==1)
                    Tlink=t1;
                else
                    Tlink=cat(2,Tlink,t1);
                end;
            end;
            if(opp1==1)
                t=Tlink;
            else
                t=cat(1,t,Tlink);
            end;
         end;
         %optimize t
         for opp11=1:degree:(kh-1)*degree+1
            for opp22=1:degree:(kw-1)*degree+1
                tt=double(t(opp11:opp11+degree-1,opp22:opp22+degree-1));
                dark_II=double(dark_I(opp11:opp11+degree-1,opp22:opp22+degree-1));
                I_rowR=double(I_row(opp11:opp11+degree-1,opp22:opp22+degree-1));
                T=tt-0.001*(dark_II-I_rowR).*tt;
                
                if(opp22==1)
                    T1=T;
                else
                    T1=cat(2,T1,T);
                end;
            end;
            if(opp11==1)
                ttt=T1;
            else
                ttt=cat(1,ttt,T1);
            end;
         end;

          LastT=max(ttt,t0);
          [h1,w1]=size(LastT);

          new1=img(:,:,1);
          new2=img(:,:,2);
          new3=img(:,:,3);
          I_1(:,:,1)=double(new1(1:h1,1:w1));
          I_1(:,:,2)=double(new2(1:h1,1:w1));
          I_1(:,:,3)=double(new3(1:h1,1:w1));


          J1(:,:,1) = uint8((I_1(:,:,1) - (1-LastT)*A)./LastT);  
          J1(:,:,2) = uint8((I_1(:,:,2) - (1-LastT)*A)./LastT);  
          J1(:,:,3) = uint8((I_1(:,:,3) - (1-LastT)*A)./LastT);

        if(lh==0&&lw==0)
            J=J1;
        end;
        if(lh==0&&lw~=0)
            right=double(img(1:h,w-lw+1:w,:));
            J=cat(2,J1,right);
            J=uint8(J);
        end;
        if(lh~=0&&lw==0)
            down=double(img(h-lh+1:h,1:w,:));
             J=cat(1,J1,down);
             J=uint8(J);
        end;    
        if(lh~=0&&lw~=0)
            right=double(img(1:h-lh,w-lw+1:w,:));
            down=double(img(h-lh+1:h,1:w,:));
            I_row=cat(2,J1,right);
            J=cat(1,I_row,down);
            J=uint8(J);
        end;  
    else
        J=img;
    end; 
    pic=J;
    img=pic;
end;
out=img;

