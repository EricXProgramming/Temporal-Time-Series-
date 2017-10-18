%%%%%%%%%%%%%%%%%%%%%%%%% ColumnNormalizer.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% This function takes the documents by topic matrix W 
% normalizes the columns of this matrix in order to make a probability
% distribution of every document over the topics. 

% Inputs: 
% W -> the documents by topics matrix.

% Output: 
% Normalized_Matrix -> the normalized documents by topics matrix.

function[Normalized_Matrix] = ColumnNormalizer(W)
sums = sum(W,1);
invSums = sums.^(-1);
invSumsMat = diag(invSums);
normedW = W*invSumsMat;
Normalized_Matrix = normedW;
end
