
for ii = 1:3
    fig = figure;
    data = zeros(500,800);

    if ii==1
        data(1,1) = 1;
    elseif ii==2
        data(1,800) = 1;
    elseif ii==3
        data(500,1) = 1;
    end
    imagesc(data);

    saveName = num2str(ii);
    exportPng4DMD(saveName, fig, 1);
    close
end