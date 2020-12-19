function img_Contours = DrawContoursAroundSegments(img, klabels)

dx = [-1, -1, 0, 1, 1, 1, 0, -1];
dy = [0, -1, -1, -1, 0, 1, 1, 1];
[m_height, m_width, m_channel] = size(img);
img_Contours = uint8(zeros(m_height, m_width));

for j = 1: m_height
    for k = 1: m_width
        np = 0;
        for i = 1: 8
            x = k+dx(1, i);
            y = j+dy(1, i);
            if (x>0&&x<=m_width&&y>0&&y<m_height)
               if (klabels(j, k)~=klabels(y, x))
                   np = np+1;
               end
            end
        end
        if (np>2)
            img_Contours(j, k) = 255;
        end
    end
end

            