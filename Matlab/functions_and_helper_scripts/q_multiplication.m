function q = q_multiplication(q1, q2)
    %% Computes the Hamilton product q = q1*q2 of two quaternions (scalar 
    % last).
    
    q(1) = q1(4)*q2(1) + q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2);
    q(2) = q1(4)*q2(2) - q1(1)*q2(3) + q1(2)*q2(4) + q1(3)*q2(1);
    q(3) = q1(4)*q2(3) + q1(1)*q2(2) - q1(2)*q2(1) + q1(3)*q2(4);
    q(4) = q1(4)*q2(4) - q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3);

end