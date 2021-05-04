function e = A2e(A)
    e = [ atan2(A(1,2,:),A(2,2,:)) ;
          -asin(A(3,2,:))          ;
          atan2(A(3,1,:),A(3,3,:)) ; ];
end