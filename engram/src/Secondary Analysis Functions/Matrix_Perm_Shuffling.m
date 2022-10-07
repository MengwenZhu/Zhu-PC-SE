function [Shuffled_Matrix]=Matrix_Perm_Shuffling(Raw_Matrix)

Shuffled_Matrix=Raw_Matrix(randperm(length(Raw_Matrix(:))));

end