A = randn(n, n); b = randn(n,1); Abk=A; pvt = 1:n;
    
% mydgetrf()
for i = 1 : n-1,
    maxind = i;
    max = abs(A(i,i));
    for t = i+1 : n
        if ( abs(A(t,i)) > max)
            maxind = t;
            max = abs(A(t,i));
        end
    end
    if (max == 0)
        disp("LU factorization failed: coefficient matrix is singular"); return;
    else
        if (maxind ~= i)
            % save pivoting information
            temps = pvt(i);
            pvt(i) = pvt(maxind);
            pvt(maxind) = temps;
            % swap rows
            tempv = A(i,:);
            A(i,:) = A(maxind,:);
            A(maxind,:) = tempv;
        end
    end
    % factorization
    for j = i+1 : n,
        A(j,i) = A(j,i)/A(i,i);
        for k = i+1 : n,
            A(j,k) = A(j,k) - A(j,i) * A(i,k);
        end
    end
end

% mydtrsm, forward substitution
y(1) = b(pvt(1))
for i = 2 : n,
    y(i) = b(pvt(i)) - sum(y(1:i-1) .* A(i, 1:i-1));
end

% mydtrsm , back substitution
x(n) = y(n) / A(n,n);
for i = n-1 : -1 : 1,
    x(i) = (y(i) - sum(x(i+1:n) .* A(i, i+1:n))) / A(i,i);
end

