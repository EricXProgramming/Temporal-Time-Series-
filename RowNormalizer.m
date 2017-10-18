%%%%%%%%%%%%%%%%%%%%%%%%%% RowNormalizer.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% This function takes the documents by topic matrix W 
% normalizes the row of this matrix in order to make a probability
% distribution of every document over the topics. 

% Inputs: 
% W -> the documents by topics matrix.

% Output: 
% Normalized_Matrix -> the normalized documents by topics matrix.

function[Normalized_Matrix] = RowNormalizer(W)
[Number_Documents,~] = size(W);  
sums = sum(W,2);
invSums = sums.^(-1);
invSumsMat = sparse(1:Number_Documents,1:Number_Documents,invSums);
Normalized_Matrix = invSumsMat*W;
end
