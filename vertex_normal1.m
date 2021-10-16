function [N, adjFList] = vertex_normal1(V)
% Vertex_normal computs per vertex unit normal 
% 
% N = vertex_normal(V)
%
% Inputs:
%   V |V| x 3 matrix of vertex positions
% 
% Outputs:
%   N a |V| x 3 matrix of unit vertex normals

Idx = knnsearch(V, V,'K', 20);
nV = size(V,1);
for index = 1: nV 
%     disp(V(index,:))
    for i= 1: 18       
         triangle_face(i,:) = [Idx(index,1),Idx(index,i+1), Idx(index,i+2)];
         s = V(triangle_face(i,:),:);
         surface_normal = surfnorm(V(triangle_face(i,:),:)');

%         vectors(i,:) = V(Idx(index,1),:) - V(Idx(index, i+1),:);
%         vector_n(i,:) = norm(vectors(i,:)); 
%         vector_norm = vector_n';
        
         face_normal(i,:) = cross(V(triangle_face(i,1),:)-V(triangle_face(i,2),:), V(triangle_face(i,1),:)-V(triangle_face(i,3),:));
%         disp(i);
         triA(i,:) = real(1/2*norm(face_normal(i,:))); 
         triArea = triA';
        
        %%%% surfnorm        
    end
    adjF{index} = [Idx(index,2), Idx(index,3), Idx(index,4), Idx(index,5), Idx(index,6), Idx(index,7), Idx(index,8)];
%     
%  
     sample_mean = sum(face_normal)/size(face_normal,1);
% 
% 
     sample_cov = ((face_normal - sample_mean) * (face_normal - sample_mean)') / size(face_normal,2);
     [M,D] = eig(sample_cov);
      eigen_trans = M';
      eigen_tran = real(eigen_trans(1:10,:));
      recon_normal = eigen_tran * (face_normal - sample_mean);
      recon_normal_mean(index,:) = sum(recon_normal) / size(recon_normal,1);
      vertex_normal(index, :) = sum(triArea * face_normal, 1)/ norm(sum(triArea * face_normal, 1));
     
     
     

    
end
% N = normalizerow(recon_normal_mean);
% N = recon_normal_mean;
N = normalizerow(vertex_normal);
adjFList = adjF';
end