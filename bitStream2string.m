function string = bitStream2string(b_hat)

% Map 0 and 1 literal values to their char representations.
char_vector = char(b_hat + '0');

% Reshape back into col-vector of binary strings.
bin_array = reshape(char_vector, 7,[]).';

% Convert binary strings into corresponding chars
% Transpose to make readable.
string = char(bin2dec(bin_array)).';