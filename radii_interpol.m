function [r] = radii_interpol(m,y1,y2,x1,x2,x)
%RADII_INTERPOL Summary of this function goes here
r = m + ((y2-y1)/(x2-x1)*x);
end

