function A = rolling_adjacency_matrix(n,n_temp)
    % Create rolling matrix 
    % Inputs:
    %   n           total nb of network nodes
    %   n_temp      max nb of nodes connected at once
    % Output:
    %   A           temporal full adjacency matrix

    % initialise A
    A=[];
    
    % create rolling matrix
    for i=1:(n-n_temp)+1
        A_temp = zeros(n,n);
        A_temp(i:i+n_temp-1,i:i+n_temp-1) = ones(n_temp,n_temp);
        A_temp(logical(eye(n)))=0;
        A = cat(3,A,A_temp);
    end

    for i=(n-n_temp)+2:n
        A_temp = zeros(n,n);
        len= n-i;
        A_temp(i:end,i:end) = ones(len+1,len+1);
        A_temp(1:(n_temp-len-1),1:(n_temp-len-1)) = ones(n_temp-len-1,n_temp-len-1);

        A_temp(1:(n_temp-len-1),i:end)=ones((n_temp-len-1),len+1);
        A_temp(i:end,1:(n_temp-len-1))=ones(len+1,(n_temp-len-1));

        A_temp(logical(eye(n)))=0;
        A = cat(3,A,A_temp);
    end
end
