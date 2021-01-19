function [eigvect,lambda] = eig_sort(C)

%function that generates eigenvalues and eigenvectors arranged in
%descending order

[eigvect_unsorted,lambda_unsorted]=eig(C);
[lambda,lambda_order]=sort(diag(lambda_unsorted),'descend');
eigvect=eigvect_unsorted(:,lambda_order);
