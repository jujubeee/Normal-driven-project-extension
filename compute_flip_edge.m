function flip_edge_energy = compute_flip_edge(have_edge,idx,PN, knn_distance, r, points_label)
    flip_edge_energy = zeros(size(idx,1), size(idx,2));
    for i = 1:size(knn_distance,1)
          for j = 1:size(knn_distance, 2)
         
              if  have_edge(i, j) == 1

                 
          
                  if (dot(PN(i,:),PN(idx(i,j),:)) < 0 && points_label(i) == points_label(idx(i,j))) || (dot(PN(i,:),PN(idx(i,j),:)) >= 0 && points_label(i) ~= points_label(idx(i,j))) 
                       flip_edge_energy(i,j) = abs(dot(PN(i,:),PN(idx(i,j),:))) * (1 - (knn_distance(i,j))^2/r^2);
                  else
                       flip_edge_energy(i,j) = 0;
                  end
             
              else
                  flip_edge_energy(i,j) = Inf;
              end
        

          end
    end
  
end
