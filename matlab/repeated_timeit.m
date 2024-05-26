function btime = repeated_timeit(func, n)
    btimes = zeros(n, 1);
    for i=1:n
        btimes(i) = timeit(func);
    end
    btime = median(btimes);
end