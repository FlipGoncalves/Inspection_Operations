function [M] = rotz(a)
% Returns rotation matrix for 3D rotation around z axis
% a - angle in radians
M = [cos(a) -sin(a) 0 0
     sin(a) cos(a) 0 0
     0 0 1 0
     0 0 0 1];
end