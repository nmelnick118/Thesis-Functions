function perimeter = polyperim(X, Y)
    L = length(X);
    perimeter = 0;
    for i = 1:L-1
        perimeter = perimeter + sqrt((X(i) - X(i+1))^2 + (Y(i)-Y(i+1))^2);
    end
    perimeter = perimeter + sqrt((X(1) - X(end))^2 + (Y(1)-Y(end))^2);
end