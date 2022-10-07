function [PV_Matrix]=PV_Correlation(Matrix1, Matrix2)

Matrix1(isnan(Matrix1))=0;
Matrix2(isnan(Matrix2))=0;

LinearizedMatrix1=nan(1,size(Matrix1,1)*size(Matrix1,2),size(Matrix1,3));
for s=1:size(Matrix1,3)
    LinearizedMatrix1(1,:,s)=reshape(Matrix1(:,:,s),1,size(Matrix1,1)*size(Matrix1,2));
end
    
LinearizedMatrix2=nan(1,size(Matrix2,1)*size(Matrix2,2),size(Matrix2,3));
for s=1:size(Matrix2,3)
    LinearizedMatrix2(1,:,s)=reshape(Matrix2(:,:,s),1,size(Matrix2,1)*size(Matrix2,2));
end


PV_Matrix=nan(1,size(Matrix1,1)*size(Matrix1,2));
for s=1:size(LinearizedMatrix1,2)
    PV_Matrix(1,s)=XCorrCalc(LinearizedMatrix1(1,s,:),LinearizedMatrix2(1,s,:));
end

end

