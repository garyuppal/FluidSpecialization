function [bvx, bvy] = brankine(rotation,radius,xp,yp)

PI = 3.141592653589793238462643383279502884197169399;

dist = sqrt(xp.^2 + yp.^2);
theta = atan2(xp,yp);

inside = (dist <= radius);
outside = (dist > radius);
factor = inside.*((rotation*dist)./(2*PI*(radius^2))) ...
       + outside.*(rotation./(2*PI*dist));

bvy = -sin(theta).*factor;
bvx = cos(theta).*factor;

end

