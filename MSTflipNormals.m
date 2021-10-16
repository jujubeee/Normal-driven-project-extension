function PN_updated = MSTflipNormals(V_ori,PN, points_label)
assert(size(PN,2) == 3)
% default parameters
u = PN(:,1);
v = PN(:,2);
w = PN(:,3);
sensorCenter = [0,0,0]; 
for k = 1: size(V_ori,1)
% for k = 1 : numel(x)
%    p1 = sensorCenter - [V_ori(k,1),V_ori(k,2),V_ori(k,3)];
%   
%    p2 = [u(k),v(k),w(k)];
   % Flip the normal vector if it is not pointing towards the sensor.
%    angle = atan2(norm(cross(p1,p2)),p1*p2');
%    if angle > pi/2 || angle < -pi/2
   if points_label(k)  == -1
       u(k) = -u(k);
       v(k) = -v(k);
       w(k) = -w(k);
       PN(k,1) = u(k);
       PN(k,2) = v(k);
       PN(k,3) = w(k);
   end
end
PN_updated = PN;
end