function optionPrices = readCSVData(fileName1, fileName2,target)
    table1=readtable(fileName1);
    table2=readtable(fileName2);
    optionPrices=[table1(table1.Var1 == target, :);
                  table2(table2.Var1 == target, :)];
    optionPrices=optionPrices(:,2:end);
    optionPrices=table2array(optionPrices);
end