function maxrelerr = discrepancy(user_data, truth_data)

if nargin ~= 2
  error('DISCREPANCY:  must have two input arrays')
end
if any(any(isnan(truth_data)))
  error('DISCREPANCY: user_data contains NaN ')
end
if any(size(user_data) ~= size(truth_data))
  error('DISCREPANCY:  two input arrays must be the same size')
end

maxrelerr = max(abs(truth_data(:)-user_data(:))./abs(truth_data(:)));

if maxrelerr > sqrt(eps)
    display('  ');
    display('Outputs do not match - UQ Challenge Code does not appear to be giving correct results ');
    display('  ');
else
    display('  ');
    display('Outputs match - UQ Challenge Code has been validated');
    display('  ');
end
