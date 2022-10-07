function [CrossCorrValue]=XCorrCalc(matrix1,matrix2)


CrossCorrValue=corr(matrix1(:),matrix2(:));


end