%% q2A - Calculate cosine matrix from q
function A = q2A(q)
    qvec = reshape(q(1:3),3,1);
    skewq = [0 -qvec(3) qvec(2); qvec(3) 0 -qvec(1) ; -qvec(2) qvec(1) 0];
    A = (q(4)^2 - sum(qvec.^2))*eye(3) + 2*qvec*(qvec') - 2*q(4)*skewq ;
end