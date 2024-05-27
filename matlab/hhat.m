function yhat = hhat(xhat, ~, ~, ~)
    xhatd = xhat(1:2);
    yd = h(xhatd);
    yhat = yd;
end
