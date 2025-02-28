function Omega = H_clock(R,O, freq)
    % Omega = A*freq + B
    % Inputs:
    %   freq                fixed hardware oscillator ticking frequency
    %   R                   hardware oscillator rate
    %   O                   hardware oscillator offset           
    % Output:
    %   Omega               virtual time  

    Omega = R.*(1./freq) + O;
    %Omega = freq + O;
end