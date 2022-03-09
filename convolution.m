function [Output] = convolution(A,Kernel)

M = size(Kernel,1)-1;

for i = 1:(size(A,1)-M)
    for j =1:(size(A,2)-M)
        Temp = A(i:i+M,j:j+M).*Kernel;
        Output(i,j)=sum(Temp(:));
    end
end
end