OM = sym('OM', [4 4]);
L = sym('L', [4 3]);
W = sym('W', [3 3]);

GA = sym('GA', [3 4]);
M = sym('M', [2 4]);

A = [OM, L, zeros(4, 8);
     zeros(3, 4), W, zeros(3,8);
     zeros(8,15)] + eye(15);
 
C = [zeros(3, 4) eye(3), eye(3), zeros(3, 5);
     GA, zeros(3,6), eye(3), zeros(3,2);
     M, zeros(2,9), eye(2)];
 
CA_mat = C;
Obsv = CA_mat;

for i=1:14
    CA_mat = CA_mat*A;
    Obsv = [Obsv; CA_mat];
end

% OM = ones(4);
% L = -1*ones(4, 3);
% W = 3*ones(3);
% 
% GA = 2*ones(3, 4);
% M = 7*ones(2, 4);
% 
% A = [OM, L, zeros(4, 8);
%      zeros(3, 4), W, zeros(3,8);
%      zeros(8,15)] + eye(15);
%  
% C = [zeros(3, 4) eye(3), eye(3), zeros(3, 5);
%      GA, zeros(3,6), eye(3), zeros(3,2);
%      M, zeros(2,9), eye(2)];
%  
% CA_mat = C;
% Obsv = CA_mat;
% 
% for i=1:14
%     disp(strcat('Calculating CA^', num2str(i)));
%     CA_mat = CA_mat*A;
%     Obsv = [Obsv; CA_mat];
% end