function [isHappy, density_a, density_b, m_happy, first5Sequences, sequenceLengths] = findHappySequence(input_sequence, n)
    % Define the transition matrix
    T = [1, 1; n, n-1];
    % Find the eigenvalues and eigenvectors
    [V, D] = eig(T);
    % Find the dominant eigenvalue and corresponding eigenvector
    [dominant_lambda, idx] = max(max(D));
    dominant_vector = V(:, idx);
    % Normalize the dominant eigenvector
    dominant_vector = dominant_vector / sum(dominant_vector);
    % Calculate the asymptotic densities
    density_a = dominant_vector(1);
    density_b = dominant_vector(2);
    % Initialize variables for happy sequence check
    sequence = input_sequence; % Start with the input sequence
    isHappy = false; % Initialize as false, will set to true if happy sequence found
    m_happy = 0; % Initialize the number of times it took for the sequence to be happy
    first5Sequences = cell(1, 5); % Initialize cell array to store first 5 new sequences
    sequenceLengths = zeros(1, 1000); % Initialize array to store the lengths of new sequences
    density_a_values = zeros(1, 10);
    density_b_values = zeros(1, 10);
    

    % Adjust m to adjust the number of steps before program terminates. 
    for m = 1:12 % A large number to represent the infinite scenario
        % Apply the mapping phi to the sequence
        new_sequence = '';
        for j = 1:length(sequence)
            if sequence(j) == '1'
                new_sequence = strcat(new_sequence, ['1', repmat('0', 1, n)]); % 1 maps to 1 followed by n zeros
            else
                new_sequence = strcat(new_sequence, ['1', repmat('0', 1, n-1)]); % 0 maps to 1 followed by n-1 zeros
            end
        end
        % Store the length of the new sequence
        sequenceLengths(m) = length(new_sequence);
        % Print the new sequence and its length
        fprintf('Step %d: New sequence = %s, Length = %d\n', m, new_sequence, sequenceLengths(m));
        
        % Calculate densities for the first 10 steps
        if m <= 10
            num_a = sum(new_sequence == '1');
            num_b = sum(new_sequence == '0');
            total_length = length(new_sequence);
            density_a_values(m) = num_a / total_length;
            density_b_values(m) = num_b / total_length;
        end

        % Store the new sequence if it's one of the first 5 iterations
        if m <= 5
            first5Sequences{m} = new_sequence;
        end
        % Check if the original sequence is a subsequence of the new sequence
        idx = strfind(new_sequence, input_sequence);
        if ~isempty(idx)
            isHappy = true;
            m_happy = m;
            % Insert brackets around the original sequence in the new sequence
            for i = length(idx):-1:1
                new_sequence = insertAfter(new_sequence, idx(i) + length(input_sequence) - 1, ']');
                new_sequence = insertBefore(new_sequence, idx(i), '[');
            end
            fprintf('Original sequence found as a subsequence: %s\n', new_sequence);
            break;
        end
        sequence = new_sequence; % Update the sequence for the next iteration
    end

    % Plotting the densities
    figure;
    plot(1:10, density_a_values, 'r-o', 1:10, density_b_values, 'b-o');
    xlabel('Step');
    ylabel('Density');
    legend('Density of 1', 'Density of 0');
    title('Asymptotic Density for the First 10 Steps');

    % If the loop completes without finding a happy sequence
    if m_happy == 0
        isHappy = false;
    end
    % Trim the sequenceLengths array to remove unused entries
    sequenceLengths = sequenceLengths(1:m_happy);
end
