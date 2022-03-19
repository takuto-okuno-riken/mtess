function b = url2cacheString(url)
    u = split(url,'?');
    u = strrep(u{1},'http://','');
    u = strrep(u,'https://','');
    u = split(u,'/');
    b = strrep(u{1},'.','_');
    for j=2:length(u)
        b = [b '-' u{j}];
    end
end
