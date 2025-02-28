function clock_param = clock_param(n)
    % Initialse node values using given distribution
    % Inputs:
    %   n                   nb of nodes
    % Output:
    %   R                   hardware oscillator rate
    %   O                   hardware oscillator offset  
    %   freq                fixed hardware oscillator ticking frequency
    %   rate                initial virtual rate
    %   offset              initial virtual offset

    R = random('Normal', 1, 0.2, [n,1]);
    R(R<0.1) = 1;
    O = random('Normal', 0, 0.2, [n,1]);
    freq = random('Normal', 1, 0.2, [n,1]);
    freq(freq<0.1) = 1;
    if sum(abs(O)>freq)>0
        O(abs(O)>freq) = 0.2*freq(abs(O)>freq);
    end

    rate = random('Normal', 1, 0.2, [n,1]);
    rate(freq<0.1) = 1;
    offset = random('Normal', 0, 0.2, [n,1]);

    dH_test = H_clock(R,O, freq);

    while sum(dH_test<0)>0
        R(dH_test<0) = random('Normal', 1, 0.2, 1);
        O(dH_test<0) = random('Normal', 0, 0.2, 1);
        freq(dH_test<0) = random('Normal', 1, 0.2, 1);

        dH_test = H_clock(R,O, freq);
    end

    %clock_param = [R, O, freq, rate, offset]
    clock_param = [R, O, freq, rate];

    
   
end