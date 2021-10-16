% function data = precomputation(V,F,lambda,type)
function data = precomputation(V,lambda,adjList, N_ori)
%     if nargin < 3
%         lambda = 0;
%         type = 'vertex';
%     elseif nargin < 4
%         type = 'vertex';
%     end
%     if (strcmp(type, 'vertex'))
    data = precomputation_vertex(V,N_ori,lambda,adjList);
%     end
    
%     if (strcmp(type, 'face'))
%         data = precomputation_face(V,F,lambda);
%     elseif (strcmp(type, 'vertex'))
%      data = precomputation_vertex(V,F,lambda);
%     end
end

function data = precomputation_vertex(V,N_ori,lambda,adjList)
% function data = precomputation_vertex(V,F,lambda)
    data.V = V;
    nV = size(V,1);
%     data.F = F;
%     data.N = vertex_normal(V); % input face normal
    data.N = N_ori;
    data.preF = []; % prefactorization of Q
    data.A = 1;
%     data.A = full(diag(massmatrix(V,F))); % face area
%     data.totalArea = sum(data.A);
    data.totalArea = nV*data.A;
%     data.cotanW = cotangent(V,F); % cotangent weight
    data.cotanW = 1;
    data.lambda = lambda;
%     adjFList = vertexFaceAdjacencyList(F);

    
    nV = size(V,1);
    
    %% construct neighboring info (Nk)
    for kk = 1:nV
        % Nk
        Nk = adjList{kk};
        
%         E_fk = [F(Nk,1) F(Nk,2); ...
%               F(Nk,2) F(Nk,3); ...
%               F(Nk,3) F(Nk,1)];
        % this is the neighbor edge of vertex
%         disp(NK);
%         disp(adjList{kk}(2));
        %%%%%% change the number of vertex 
        E_fk = [kk, Nk(1);
            kk, Nk(2);
            kk, Nk(3);
            kk, Nk(4);
            kk, Nk(5);
%             Nk(1), Nk(2);
%             Nk(1), Nk(3);
%             Nk(2), Nk(3);
             kk, Nk(6);
             kk, Nk(7);
%             kk, Nk(8);
%              kk, Nk(9);
%             kk, Nk(10)
%             Nk(1), Nk(2);
%             Nk(2), Nk(3); 
%             Nk(2), Nk(4);
%             Nk(3), Nk(4);
%             Nk(4), Nk(5);
%             Nk(1), Nk(3);
%             Nk(1), Nk(4);
%             Nk(1), Nk(5);
                ];
       %%%% should be vector; 
       W_fk = ones(7,1);

          
        % save info
        data.E_f{kk} = E_fk;
        data.W_f{kk} = diag(W_fk);
        data.dV_f{kk} = (V(E_fk(:,2),:) - V(E_fk(:,1),:))';
    end
    
    %% precomputation for the global step (Q1 and K1)
    QIJV = zeros(nV*3*4*4,3); % construct a long enough list for Q1
    KIJV = zeros(nV*18*3*4,3); % construct a long enough list for K1
    QIdx = 1;
    KIdx = 1;
    for kk = 1:nV
        
        E_fk = data.E_f{kk};
        W_fk = diag(data.W_f{kk});
        
        nE = size(E_fk,1);
        
        Qi = [E_fk(:,1); E_fk(:,2); E_fk(:,1); E_fk(:,2)];
        Qj = [E_fk(:,2); E_fk(:,1); E_fk(:,1); E_fk(:,2)];
         Qv = [W_fk; W_fk; -W_fk; -W_fk];  
%         Qv = ones(size(Qi,1), size(Qi,2));
        QIJV(QIdx:QIdx+4*nE-1,:) = [Qi,Qj,Qv];
        QIdx = QIdx + 4*nE;
        
        Ki = [repmat(1,nE*2,1); repmat(2,nE*2,1); repmat(3,nE*2,1)];
        Kj = repmat([E_fk(:,1); E_fk(:,2)], 3, 1);
%         Kv = [1 .* V(E_fk(:,1),:)-1 .* V(E_fk(:,2),:); 1 .* V(E_fk(:,2),:)- 1 .* V(E_fk(:,1),:)];
        Kv = [W_fk .* V(E_fk(:,1),:)-W_fk .* V(E_fk(:,2),:); W_fk .* V(E_fk(:,2),:)-W_fk .* V(E_fk(:,1),:)];
        Kv = Kv(:);
        KIJV(KIdx:KIdx+6*nE-1, :) = [Ki+9*(kk-1), Kj, Kv];
        KIdx = KIdx + 6*nE;
    end
    QIJV(QIdx:end,:) = [];
    KIJV(KIdx:end,:) = [];
    
    data.Q1 = sparse(QIJV(:,1),QIJV(:,2),QIJV(:,3),size(V,1),size(V,1));
    
    data.K1 = sparse([KIJV(:,1); KIJV(:,1)+3; KIJV(:,1)+6], ...
                [KIJV(:,2); KIJV(:,2)+nV; KIJV(:,2)+nV+nV], ...
                [KIJV(:,3);KIJV(:,3);KIJV(:,3)], ...
                9*nV,3*nV);  
            
    data.LHS = data.Q1 /2;

end

% function data = precomputation_face(V,F,lambda)
%     data.V = V;
%     data.F = F;
%     data.N = faceNormals(V,F); % input face normal
%     data.preF = []; % prefactorization of Q
%     data.A = doublearea(V,F) / 2; % face area
%     data.totalArea = sum(data.A);
%     data.cotanW = cotangent(V,F); % cotangent weight
%     data.lambda = lambda;
%     
%     nF = size(F,1);
%     
%     %% construct neighboring info (Nk)
%     for kk = 1:size(F,1)
%         % Nk
%         Nk = kk;
%         
%         E_fk = [F(Nk,1) F(Nk,2); ...
%               F(Nk,2) F(Nk,3); ...
%               F(Nk,3) F(Nk,1)];
%         
%         W_fk = [data.cotanW(Nk, 3); ...
%               data.cotanW(Nk, 1); ...
%               data.cotanW(Nk, 2)];
%           
%         % save info
%         data.E_f{kk} = E_fk;
%         data.W_f{kk} = diag(W_fk);
%         data.dV_f{kk} = (V(E_fk(:,2),:) - V(E_fk(:,1),:))';
%     end
%     
%     %% precomputation for the global step (Q1 and K1)
%     QIJV = zeros(size(F,1)*3*4,3); % construct a long enough list for Q1
%     KIJV = zeros(size(F,1)*18,3); % construct a long enough list for K1
%     QIdx = 1;
%     KIdx = 1;
%     for kk = 1:size(F,1)
%         E_fk = data.E_f{kk};
%         W_fk = diag(data.W_f{kk});
%         
%         nE = size(E_fk,1);
%         
%         Qi = [E_fk(:,1); E_fk(:,2); E_fk(:,1); E_fk(:,2)];
%         Qj = [E_fk(:,2); E_fk(:,1); E_fk(:,1); E_fk(:,2)];
%         Qv = [W_fk; W_fk; -W_fk; -W_fk];        
%         QIJV(QIdx:QIdx+4*nE-1,:) = [Qi,Qj,Qv];
%         QIdx = QIdx + 4*nE;
%         
%         Ki = [repmat(1,nE*2,1); repmat(2,nE*2,1); repmat(3,nE*2,1)];
%         Kj = repmat([E_fk(:,1); E_fk(:,2)], 3, 1);
%         Kv = [W_fk .* V(E_fk(:,1),:)-W_fk .* V(E_fk(:,2),:); W_fk .* V(E_fk(:,2),:)-W_fk .* V(E_fk(:,1),:)];
%         Kv = Kv(:);
%         KIJV(KIdx:KIdx+6*nE-1, :) = [Ki+9*(kk-1), Kj, Kv];
%         KIdx = KIdx + 6*nE;
%     end
%     fprintf("QIdx %d, QIJV length: %d\n", QIdx, length(QIJV));
%     fprintf("KIdx %d, KIJV length: %d\n", KIdx, length(KIJV));
%     QIJV(QIdx:end,:) = [];
%     KIJV(KIdx:end,:) = [];
%     
%     data.Q1 = sparse(QIJV(:,1),QIJV(:,2),QIJV(:,3),size(V,1),size(V,1));
%     
%     data.K1 = sparse([KIJV(:,1); KIJV(:,1)+3; KIJV(:,1)+6], ...
%                 [KIJV(:,2); KIJV(:,2)+nF; KIJV(:,2)+nF+nF], ...
%                 [KIJV(:,3);KIJV(:,3);KIJV(:,3)], ...
%                 9*nF,3*nF);  
%             
%     data.LHS = data.Q1 /2; 
% end
            
    
