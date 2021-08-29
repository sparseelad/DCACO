function X=vec_to_smatrix(z,n)


%====================================================
%This Script constructs the symmetric matrix W with 
%diagonal elements being all 1 from a given vector w.


%====================================================



if length(z)~=(n^2-n)/2
    
    error('Dimension of the input vector is wrong');

else
    
    X=eye(n);
    
    for i=1:(n-1)
        for j=(i+1):n
            s=j-i+((2*n-i)*(i-1))/2;
            X(i,j)=z(s);
            X(j,i)=z(s);
            
        end
    end
    
end

end

