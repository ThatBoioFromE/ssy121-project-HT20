function constellation_vec = generateConstellation(type)
% All constellations normalized to max distance = 1!

switch(type)
    case 'BPSK'
        constellation_vec = [1 -1]; % BPSK
    case '4-QAM'
        constellation_vec = [(1 + 1i), (1 - 1i), (-1 -1i), (-1 + 1i)]/sqrt(2);% 4-QAM
    case '8-PSK'
        constellation_vec = [sqrt(2), (1 + 1i), sqrt(2)*1i, (-1 + 1i), -sqrt(2),...
            (-1 -1i), -sqrt(2)*1i, (1 - 1i)]/sqrt(2);  % 8-PSK
    case '16-QAM'
        % This is tricky because, although the points lie on three
        % concentric circles, they are not evenly spaced -> assign using
        % some kind of grid rather than rotation.
        
        % Generate matrix of 16-QAM points.
        order = 16  % Number of QAM constellation points.
        d_min = 2;  % Minimum distance between symbols.
        constx = 0:d_min:(sqrt(order)-1)*d_min;
        constx = repmat(constx, sqrt(order),1);
        consty = 0:d_min:(sqrt(order)-1)*d_min;
        consty = repmat(consty', 1, sqrt(order));
        const = constx+1i*consty; % Make each point into complex number.
        const = const - (3+3i); % Translate matrix to center on origin.
        % Next step: organise constellation matrix elements into vector sorted in
        % gray sequence.

        constellation_vec = reshape(const, [1 16])./3; % Normalize, giving max magnitude = 1.
        % plot(const, 'ob')
        % grid on
        
    otherwise
        disp('Unrecognized constellation type. Provide a valid string')
end