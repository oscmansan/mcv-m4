function [angle] = get_angle_between_two_lines(l, m)
    l = [l(1)/l(3);l(2)/l(3);1];
    m = [m(1)/m(3);m(2)/m(3);1];
    omega = [1, 0, 0;
             0, 1, 0;
             0, 0, 0];
    angle = acos(dot((omega*l),m)/(sqrt(dot(omega*l,l))*sqrt(dot(omega*m,m))));
end