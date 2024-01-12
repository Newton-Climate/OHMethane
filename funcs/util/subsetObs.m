function out = subsetObs(in, inds)

out = struct();
fields = fieldnames(in);
for i=1:length(fields)
field = char(fields(i));
out.(field) = in.(field)(inds);
end
