function PN = filpNormals_cheated(V_ori,PN, N_ori)
% assert(size(PN,2) == 3)
% % default parameters
% u = PN(:,1);
% v = PN(:,2);
% w = PN(:,3);
% sensorCenter = [0,0,0]; 
PN_flipped = - PN;
 for k = 1: size(V_ori,1)
% % for k = 1 : numel(x)
%    p1 = sensorCenter - [V_ori(k,1),V_ori(k,2),V_ori(k,3)];


minvalue= min(sqrt(sum((N_ori(k,:) - PN(k,:)).^2,2)), sqrt(sum((N_ori(k,:) - PN_flipped(k,:)).^2,2)));
if minvalue == sqrt(sum((N_ori(k,:) - PN_flipped(k,:)).^2,2))
    PN(k, :) = -PN(k,:);
end
%    p2 = [u(k),v(k),w(k)];
%    % Flip the normal vector if it is not pointing towards the sensor.
%    angle = atan2(norm(cross(p1,p2)),p1*p2');
% %    if angle > pi/2 || angle < -pi/2
%    if angle >= -pi/2 && angle <= pi/2
%        u(k) = -u(k);
%        v(k) = -v(k);
%        w(k) = -w(k);
%        PN(k,1) = u(k);
%        PN(k,2) = v(k);
%        PN(k,3) = w(k);
%    end
end
PN = PN;
end