x = gen_lat_sq(600, 6);
x = [x, fliplr(x)];

for r = 1:size(x,1)
    % ref: https://www.mathworks.com/matlabcentral/answers/21-how-do-i-convert-a-numerical-vector-into-a-comma-delimited-string
    allOneString = sprintf('%.0f,' , x(r,:));
    allOneString = ['[', allOneString(1:end-1), ']'];
    
    c{r,1} = allOneString;
    c{r,2} = 0;
end

%Create CSV
tab = cell2table(c, 'VariableNames', {'cond_seq', 'claimed'});
writetable(tab,'latsq.csv');