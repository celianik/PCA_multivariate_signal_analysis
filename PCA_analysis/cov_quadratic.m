% create a matrix
S = [1 2; 2 9];

% weights along each dimension
wi = -2:.1:2;

quadform = zeros(length(wi));

for i = 1:length(wi)
  for j = 1:length(wi)
    
    % column vector
    w = [wi(i) wi(j)]';
    norm_factor = w'*w;
    quadform(i,j) = w'*S*w/norm_factor;
    
  end
end

figure(1), clf

surf(wi,wi,quadform'), shading interp
title('Visual represantation of quadratic form of S')
xlabel('W1'), ylabel('W2'), zlabel('energy')
rotate3d on

[eig_vecs, eig_val] = eig(S);

hold on
plot3([0 eig_vecs(1,1)],[0 eig_vecs(2,1)], [1 1]*max(quadform(:)),'r','linew',3)
plot3([0 eig_vecs(1,2)],[0 eig_vecs(2,2)], [1 1]*max(quadform(:)),'r','linew',3)