%% A2q - Calculate quaternion (scalar last) from cosine matrix via Stanley's method
function q = A2q(A)
    trA = trace(A);
    beta4 = 0.25*(1+trA);
    beta_vec = 0.25*(1 + 2*diag(A) - trA);
    a = [1; 2; 3];
    if all(beta4 > beta_vec) % if beta4 is greatest
        beta4 = sqrt(beta4);
        is = mod(a,3)+1;
        js = mod(a+1,3)+1;
        beta_vec = 0.25/beta4*(A(sub2ind([3 3],is,js)) - A(sub2ind([3 3],js,is)));
    else
        [bmax, maxi] = max(beta_vec);
        bmax = sqrt(bmax);
        i = mod(maxi,3)+1;
        j = mod(maxi+1,3)+1;
        beta4 = 0.25/bmax *(A(i,j) - A(j,i));
        is = a(a ~= maxi);
        beta_vec(is) = 0.25/bmax * (A(maxi,is)' + A(is,maxi));
    end
    q = [beta_vec; beta4];
end