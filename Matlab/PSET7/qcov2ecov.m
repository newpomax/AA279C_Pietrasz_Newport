function ecov = qcov2ecov(qcov, q)
    G = zeros(3,4);
    G(3,1) = -(q(3)+q(2))/((q(3)+q(2))^2 + (q(4)+q(1))^2) + (q(3)-q(2))/((q(3)-q(2))^2 + (q(4)-q(1))^2);
    G(3,2) = (q(4)+q(1))/((q(3)+q(2))^2 + (q(4)+q(1))^2) - (q(4)-q(1))/((q(3)-q(2))^2 + (q(4)-q(1))^2);
    G(3,3) = (q(4)+q(1))/((q(3)+q(2))^2 + (q(4)+q(1))^2) + (q(4)-q(1))/((q(3)-q(2))^2 + (q(4)-q(1))^2);
    G(3,4) = -(q(3)+q(2))/((q(3)+q(2))^2 + (q(4)+q(1))^2) - (q(3)-q(2))/((q(3)-q(2))^2 + (q(4)-q(1))^2);
    meq = 2/sqrt(1-4*(q(2)*q(3)+q(1)*q(4))^2);
    G(2,1:4) = flip(q)*meq;
    G(1,1:4) = flip(G(3,1:4));
    ecov = G*qcov*(G.');
end