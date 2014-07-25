function array = EigSort(eigs)

for i=length(eigs):-1:1
    for j=1:i-1
        if abs(eigs(j)) < abs(eigs(j+1)) || ...
           abs(eigs(j))==abs(eigs(j+1)) &&  real(eigs(j)) < real(eigs(j+1)) || ...
           abs(eigs(j))==abs(eigs(j+1)) &&  real(eigs(j)) == real(eigs(j+1)) && ...
           imag(eigs(j)) < imag(eigs(j+1))
            temp=eigs(j+1);
            eigs(j+1)=eigs(j);
            eigs(j)=temp;
        end        
    end
end

array=eigs;

end