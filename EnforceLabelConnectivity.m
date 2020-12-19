function nlabels = EnforceLabelConnectivity(img_Lab, labels, K)

dx = [-1, 0, 1, 0];
dy = [0, -1, 0, 1];
[m_height, m_width, m_channel] = size(img_Lab);
[M, N] = size(labels);
SUPSZ = (m_height*m_width)/K;  
nlabels = (-1)*ones(M, N);

label = 1;
adjlabel = 1;
xvec = zeros(m_height*m_width, 1);
yvec = zeros(m_height*m_width, 1);
m = 1;
n = 1;

for j = 1: m_height
    for k = 1: m_width
        if (0>nlabels(m, n))
            nlabels(m, n) = label;
            xvec(1, 1) = k;
            yvec(1, 1) = j;
            for i = 1: 4
                x = xvec(1, 1)+dx(1, i);
                y = yvec(1, 1)+dy(1, i);
                if (x>0 && x<=m_width && y>0 && y<=m_height)
                    if (nlabels(y, x)>0)
                        adjlabel = nlabels(y, x);
                    end
                end
            end
            
            count = 2;
            c = 1;
            while (c<=count)
                for i = 1: 4
                    x = xvec(c, 1)+dx(1, i);
                    y = yvec(c, 1)+dy(1, i);
                    if (x>0 && x<=m_width && y>0 && y<=m_height)
                        if (0>nlabels(y, x) && labels(m, n)==labels(y, x))
                            xvec(count, 1) = x;
                            yvec(count, 1) = y;
                            nlabels(y, x) = label;
                            count = count+1;
                        end
                    end
                end
                c = c+1;
            end
            
            if (count<(SUPSZ/4))
                for c = 1: (count-1)
                    nlabels(yvec(c, 1), xvec(c, 1)) = adjlabel;
                end
                
                label = label-1;
            end
            
            label = label+1;
            
        end
       
        n = n+1;
        if (n>m_width)
            n = 1;
            m = m+1;
        end
    
    end
end

end


