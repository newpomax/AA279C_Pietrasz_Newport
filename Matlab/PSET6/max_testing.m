w = ones(7,1);
A = [0.9945    0.1044    0.0000 ;
   -0.1044    0.9945   -0.0000 ;
   -0.0000    0.0000    1.0000 ] ;
M = [[-0.1; 0.8; 0]/sqrt(0.064+0.01), A.', A.'];
V = M;
qmethod(w,M,V)
function q_est = qmethod(w,M,V)
    w = reshape(w,1,[]); % row vec
    W = M.*(ones(size(M,1),1)*sqrt(w)); % Wi = sqrt(wi)*Mi
    U = V.*(ones(size(V,1),1)*sqrt(w)); % Ui = sqrt(wi)*Vi
    B = W*(U.');
    S = B + B.' ;
    Z = [ B(2,3)-B(3,2) ; B(3,1)-B(1,3) ; B(1,2)-B(2,1) ];
    sig = trace(B);
    K = [ S-sig*eye(size(S)) ; Z.' ];
    K = [K, [Z; sig]];
    [eVec, e] = eig(K);
    e = diag(e); % take to vector
    [~, maxi] = max(e);
    q_est = real(eVec(:,maxi));
    q_est = q_est/norm(q_est); % normalize
end