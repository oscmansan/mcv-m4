function [line] = get_line_from_points(x, y)
    coeff = polyfit([x(1), y(1)], [x(2), y(2)], 1); 
    line = [coeff(1) -1 coeff(2)];
end