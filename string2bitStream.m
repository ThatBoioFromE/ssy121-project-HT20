function b = string2bitStream(string)

% string = 'supercalifragilisticexpialidocious';

% Convert string chars into col-vector of binary strings.
bin_array = dec2bin(string, 7);

% Transpose array (due to the way matrix elements are enumerated, this will let reshape work).
% Subtract '0': ie. map it from ASCII representation of '0' to actual 0 (ie. NULL in ASCII). This automatically convert it to a double-array.
% Reshape into row-vector. 
% '[]' lets 'reshape' automatically calculate row size.
b = reshape(bin_array.'-'0', 1, []);

