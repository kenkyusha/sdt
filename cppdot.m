function r = cppdot( x, y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if (length(x) ~= length(y))
        return; 
    else
        r = 0;
        for i = 1:length(x)
            r = r + x(i)*y(i);
        end
    end


end

