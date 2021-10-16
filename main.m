clc; clear all; close all;

%% load input mesh
% [V,F] = readOBJ('spot.obj');
[V_ori,F_ori] = readOBJ('spot.obj');
V = readOBJ('spot.obj');
nV = size(V,1); % number of vertices


%%%%%%%%%%%%%%%%%%%%
%  x = V(:,1) ; 
%  y = V(:,2) ; 
%  z = V(:,3) ; 
%  dt = delaunayTriangulation(x,y, z) ;
%  tri = dt.ConnectivityList ;
%  x = dt.Points(:,1) ;
%  y = dt.Points(:,2) ; 
%  trisurf(tri,x,y,z)
%  shading interp
%  view(2)

%%%%%%%%%%%%%%%%%%%%%
% Idx = knnsearch(V, V,'K', 20);
% for index = 1: nV 
% %     disp(V(index,:))
%     for i= 1: 18       
%         triangle_face(i,:) = [Idx(index,1),Idx(index,i+1), Idx(index,i+2)];
% %         s = V(triangle_face(i,:),:);
% %         surface_normal = surfnorm(V(triangle_face(i,:),:)');
%         normal(i,:) = cross(V(triangle_face(i,1),:)-V(triangle_face(i,2),:), V(triangle_face(i,1),:)-V(triangle_face(i,3),:));
%         %%%% surfnorm
%         
%     end
%  
%     sample_mean = sum(normal)/size(normal,1);
%     disp(sample_mean);
% 
%     sample_cov = ((normal - sample_mean) * (normal - sample_mean)') / size(normal,2);
%     [M,D] = eig(sample_cov);
%      eigen_trans = M';
%      eigen_tran = eigen_trans(1:10,:);
%      recon_normal = eigen_tran * (normal - sample_mean);
%      recon_normal_mean(index,:) = sum(recon_normal) / size(recon_normal,1);
% %      normal = [];
%     
% end
%%%%%%%%%%%%%%%%%%%%%

%% compute desired normals (t)
% compute unit normal of each vertex
% N = per_vertex_normals(V,F); 

%%%%%%% method 1%%%%%%%%
% [N, adjFList] = vertex_normal1(V_ori);
% N = flipNormals(V_ori, N);
% plotMesh(V_ori, F_ori, 'fv', (N+1)/2);


%%%%%%method 2 %%%%%%%%%
% [N, adjFList] = vertex_normal2(V_ori);
% N = flipNormals(V_ori, N);
% plotMesh(V_ori, F_ori, 'fv', (N+1)/2);

%%%%%%%method 3 PCA%%%%%%%%
% [N,adjFList] = vertex_normal(V);
% N = flipNormals(V_ori, N);
% plotMesh(V, F_ori, 'fv', (N+1)/2);


%%%%%%method 4%%%%%%%%%%

 Num_porint = 10;
 %%%%%%%MST
idx = knnsearch(V_ori,V_ori,'K',Num_porint+1);
%  pcshow(V)
%  hold on
%  scatter3(V(59,1),V(59,2),V(59,3))
%  scatter3(V(1371,1),V(1371,2),V(1371,3))

%  adjF = zeros(size(V,1), 1);
 idx(:,1) = [];
% % estimate normals
 PN = zeros(size(V_ori,1),3);
 for ii = 1:size(V_ori,1)
     Ni = idx(ii,:);
     dV = V_ori(ii,:) - V_ori(Ni,:);
     [U,sigma,~] = svd(dV'*dV);
     PN(ii,:) = U(:,end)';
     
     adjF{ii} = [idx(ii,1),idx(ii,2),idx(ii,3),idx(ii,4),idx(ii,5),idx(ii,6),idx(ii,7),idx(ii,8),idx(ii,9),idx(ii,10)];
 
 end
 r = 0.12;
[knn_distance, have_edge] = compute_have_edge(idx,V,r); 
 Label = [1,-1];
 points_label = ones(size(V,1), 1);
 %%%% random assign the first point a label


flip_edge_energy = compute_flip_edge(have_edge,idx,PN, knn_distance, r, points_label);
smallest_energy = 1000;
smallest_start = -1;
% for start_point = 1: size(V,1)
%     [sum_of_energy, points_label] = MST(start_point,flip_edge_energy,V, idx,points_label);
%     disp(sum_of_energy)
%     disp(start_point)
% if sum_of_energy < smallest_energy
%     smallest_energy = sum_of_energy;
%     smallest_start = start_point;
% end
% 
% end
[sum_of_energy, points_label] = MST(1,flip_edge_energy,V, idx,points_label,Num_porint, PN, knn_distance, r);
 adjList = adjF';
PN_updated = MSTflipNormals(V_ori, PN, points_label);

% % plotMesh(V_ori,F_ori,'fv',(PN+1)/2);
% PN = flipNormals(V_ori, PN);

 
% plotMesh(V_ori, F_ori, 'fv', (PN+1)/2);
% %%%%%%%%%%%%%%%
%%
N_ori = per_vertex_normals(V_ori,F_ori); 
% PN = filpNormals_cheated(V_ori,PN, N_ori);
% [N, adjList] = vertex_normal(V);
% plotMesh(V_ori,F_ori,'fv', (N_ori+1)/2);
err = sqrt(sum((PN_updated  - N_ori).^2,2));
for i = 1:size(err,1)
    if  err(i,1) > 1                 
        disp(i)
    end
end
plotMesh(V,F_ori, 'fv', err);
colorbar;
%%
%%%%%%%%%%%%%%%%%%

% ptCloud = pointCloud(V_ori);
% Num_point = 10;
% PN = pcnormals(ptCloud, Num_point);
% % plotMesh(V_ori,F_ori,'fv',(PN+1)/2);
% u = PN(:,1);
% v = PN(:,2);
% w = PN(:,3);
% sensorCenter = [0,0,0]; 
% for k = 1: size(V_ori,1)
% % for k = 1 : numel(x)
%    p1 = sensorCenter - [V_ori(k,1),V_ori(k,2),V_ori(k,3)];
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
% end
% plotMesh(V_ori, F_ori, 'fv', (PN+1)/2);
% figure
% pcshow(V)
% color = (N+1)/2;
% S = repmat([25],2930,1);
% s = S(:);
% figure
% % scatter3(V(:,1), V(:,2), V(:,3), [], color);
% % scatter3(V(:,1), V(:,2), V(:,3), S, color);
% pcshow(V)
% title('Estimated Normals of Point Cloud')
% hold on
% x = V(:,1);
% y = V(:,2);
% z = V(:,3);
% u = N(:,1);
% v = N(:,2);
% w = N(:,3);
% % quiver3(x,y,z,u,v,w);
% quiver3(x,y,z,N(:,1), N(:,2), N(:,3));
% sensorCenter = [0,0,0]; 
% for k = 1: size(V,1)
% % for k = 1 : numel(x)
%    p1 = sensorCenter - [x(k),y(k),z(k)];
%    p2 = [u(k),v(k),w(k)];
%    % Flip the normal vector if it is not pointing towards the sensor.
%    angle = atan2(norm(cross(p1,p2)),p1*p2');
% %    if angle > pi/2 || angle < -pi/2
%    if angle >= -pi/2 && angle <= pi/2
%        u(k) = -u(k);
%        v(k) = -v(k);
%        w(k) = -w(k);
%        N(k,1) = u(k);
%        N(k,2) = v(k);
%        N(k,3) = w(k);
%    end
% end
% 
% figure
% pcshow(V)
% title('Adjusted Normals of Point Cloud')
% hold on
% quiver3(x, y, z, u, v, w);
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calucate point normal here
% plt(V);
%% precomputation
lambda = 1;
% dats = precomputation(V,lambda);
data = precomputation(V,lambda, adjList, N_ori);


% snap to closest axis
[~,cIdx] = max(abs(N_ori), [],2);
for ii = 1:size(N_ori,1)
    n = N_ori(ii,cIdx(ii));
    N_ori(ii,:) = 0;
    N_ori(ii,cIdx(ii)) = sign(n);
end
% [~,cIdx] = max(abs(N), [],2);
% for ii = 1:size(N,1)
%     n = N(ii,cIdx(ii));
%     N(ii,:) = 0;
%     N(ii,cIdx(ii)) = sign(n);
% end
% data.t = N;
data.t = N_ori;

%% start optimization
U = V; % output vertex positions

% optimization parameter
tolerance = 1e-4; 
maxIter = 100;

% we have to pin down at least one vertex to avoid the mesh flies away
% b = F(1,1); 
% bc = U(b,:);
bc = V(1,:);
b = 1;
objHis = [];
UHis = zeros(size(V,1), size(V,2), maxIter+1);
UHis(:,:,1) = U;

for iter = 1:maxIter
    
    % local step
    [RAll, objVal, data] = fitRotation_normal(U, data);
    
    % save optimization info
    objHis = [objHis objVal];
    UHis(:,:,iter+1) = U; 
    
    % global step
    Rcol = reshape(permute(RAll,[2,1,3]),1,nV*3*3);
    RHScol = data.K1' * Rcol';
    RHS = reshape(RHScol,size(RHScol,1)/3, 3);
    UPre = U;
    
    % note this min_quad_with_fixed is the same as solving Q*V = K'*R'
%     [U,data.preF] = min_quad_with_fixed(data.LHS,RHS,b,bc,[],[],data.preF);
    [U,data.preF] = min_quad_with_fixed(data.LHS,RHS,b,bc,[],[],[]);
    
    % plot
    if mod(iter-1,5) == 0
        figure(1)
        subplot(1,2,1)
%         tsurf(F,V);
        pcshow(V);
        axis equal
        subplot(1,2,2)
%         tsurf(F,U);
%         pcshow(V);
        pcshow(U);
        axis equal
        drawnow
    end

    % check whether to stop optimization
    dU = sqrt(sum((U - UPre).^2,2));
    dUV = sqrt(sum((U - V).^2,2));
    reldV = max(dU) / max(dUV);
    fprintf('iter: %d, obj: %d, reldV: %d\n', ...
        [iter, objVal, reldV]);
    
end
% writeOBJ('outupt.obj', U, F)
writeOBJ('outupt.obj', U,[]);


