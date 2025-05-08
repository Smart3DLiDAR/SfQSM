function [x, y] = bezir_n(points, dots)
    % points format:
    %   size = [2, n + 1], including starting point & end point.
    % dots:
    %   generate dots points in the bezir curve.
    
    % check points size
    n = length(points) - 1;
    if n < 1
        x = []; y = [];
        return;
    elseif n == 1
        x = points(1, :);
        y = points(2, :);
        return;
    else
        % calculate k
        k = [1 1];
        for i = 1 : (n - 1)
            k = conv(k, [1 1]);
        end
        
        % calculate t
        t = linspace(0, 1, dots);
        
        % x and y are obtained by the sum of n + 1 polynomials.
        x = zeros(size(t));
        y = zeros(size(t));
        for i = 1 : (n + 1)
            j = i - 1;
            ft = (1 - t).^(n - j).*(t.^j);
            x = x + points(1, i) * k(i) * ft;
            y = y + points(2, i) * k(i) * ft;
        end
    end
end
