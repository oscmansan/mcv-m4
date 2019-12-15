function angle = angle_between_lines(l, m)
l = [l(1)/l(3); l(2)/l(3)];
m = [m(1)/m(3); m(2)/m(3)];
angle = acos((l'*m)/(sqrt(l'*l)*sqrt(m'*m)));
end