function S = SumMatrixList(C)

S=C{1};
for i=2:length(C)
    S = S + C{i};
end

end