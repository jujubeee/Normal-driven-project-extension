function [knn_distance, have_edge] = compute_have_edge(idx,V,r)
     knn_distance = zeros(size(idx,1), size(idx,2));
     knn_dis = zeros(size(idx,1), size(idx,2));
     for i = 1:size(idx, 1)
         for j = 1:size(idx,2)
             knn_dis(i,j) = norm(V(i, :) - V(idx(i,j), :));
             if norm(V(i, :) - V(idx(i,j), :)) <= r
                knn_distance(i, j) = norm(V(i, :) - V(idx(i,j), :));
             end
         end

     end
     %%%% calculate if the two points have edge connect them 
      have_edge = zeros(size(idx,1), size(idx,2));
     for i = 1:size(knn_distance,1)
         for j = 1:size(knn_distance, 2)
             if knn_distance(i,j) ~= 0
                  have_edge(i, j) = 1;  
%                     have_edge{i}=[];
%                     disp(have_edge)
%                    have_edge{i} = [have_edge{i},idx(i,j)];
             end

         end
     end
    
end
