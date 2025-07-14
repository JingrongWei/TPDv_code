function [Ai, Bi, BTi] = mgmat(A, Aiold, Res, Pro)
    level = length(Aiold);
    
    Ai = cell(level, 1);
    Ai{level} = A;
    if level == 1
        Bi{1} = tril(Ai{1}); 
        BTi{1} = transpose(Bi{1}); 
        Res = []; Pro = [];
    end
    for j = level:-1:2
        Ai{j-1} = Res{j}*Ai{j}*Pro{j-1};           % Ac = Res*Af*Pro
        % GS smoother
        Bi{j} = tril(Ai{j});        % Forward Gauss-Seidel   B = D+L
        BTi{j} = triu(Ai{j});       % Backward Gauss-Seidel BT = D+U     
    end
end