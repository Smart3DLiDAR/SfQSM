function edge_num = plotConnectivity(pts, A, sizee,colore)
%根据邻接关系连接各个边
    A( logical(eye(size(A)))) = 0;
    tmax = max(max(A))*1.0;
    edge_num = 0;
    
    if nargin<3            
        for i=1:(size(A,1)-1)
            for j=(i+1):size(A,2)
                if( A(i,j)>0 )
                    idx = [i;j];
                    line( pts(idx,1),pts(idx,2),pts(idx,3),'LineWidth', 150*A(i,j)/tmax, 'Color', [A(i,j)/tmax,.0,.0]);
                    %a.Color(4) = 0.3;
                    edge_num = edge_num +1;
                end
            end
        end
    else
        for i=1:(size(A,1)-1)
            for j=(i+1):size(A,2)
                if( A(i,j)>0 )
                    idx = [i;j];
                    line( pts(idx,1),pts(idx,2),pts(idx,3), 'LineWidth', sizee, 'Color', colore);
                    %a.Color(4) = 0.3;
                    %colorid = randi(256,1,1);
                    %colorline = colore(colorid,:);
                    %line( pts(idx,1),pts(idx,2),pts(idx,3), 'LineWidth', sizee, 'Color', colorline);
                    edge_num = edge_num +1;
                end
            end
        end
    end
end
