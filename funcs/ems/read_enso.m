function out = read_enso(filename, sYear, eYear)

ds = csvread(filename, 1,0);
out.time = ds(:,1);
out.index = ds(:,end);
end
