function [output] = hyperrectangle_area(input)

intervals = zeros(size(input,1),1);
for i = 1:size(input,1)
intervals(i) = abs(input(i,2) - input(i,1));
end

output.area = prod(intervals);
output.intervals = intervals;
end