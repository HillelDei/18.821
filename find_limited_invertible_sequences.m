function find_limited_invertible_sequences()
    N = input('Enter the value for N for the transformation: ');
    forward_times = input('Enter the number of times the forward transformation should be done: ');
    backward_times = input('Enter the number of times the backward transformation should be done: ');

    max_trials = 3;
    sequence_length = 50;
    
    unique_sequences = cell(0);

    for trial = 1:max_trials
        % Generate initial sequence
        seq = generate_initial_sequence(N, sequence_length);
        
        original_seq = seq;

        % Apply forward transformations
        for j = 1:forward_times
            seq = apply_forward_transformation(seq, N);
        end
        
        % Apply backward transformations
        for j = 1:backward_times
            seq = apply_backward_transformation(seq, N);
            if isempty(seq)
                break;
            end
        end
        
        if isequal(seq, original_seq)
            unique_sequences{end+1} = original_seq;
            plot_sequence(original_seq);  % Plot the sequence
        end
    end

    disp('Unique sequences that regenerate themselves:');
    for i = 1:length(unique_sequences)
        disp(unique_sequences{i});
    end
end

function seq = generate_initial_sequence(N, len)
    unit1 = [1, zeros(1, N)];
    unit2 = [1, zeros(1, N-1)];
    
    seq = [];
    last_choice = 0;
    counter_unit2 = 0;
    
    while length(seq) < len
        choice = randi(2);
        
        if last_choice == 1 && choice == 1
            choice = 2;
        end
        
        if counter_unit2 >= 2
            choice = 1;
        end
        
        if choice == 1
            seq = [seq, unit1];
            counter_unit2 = 0;
        else
            seq = [seq, unit2];
            counter_unit2 = counter_unit2 + 1;
        end
        
        last_choice = choice;
    end
    
    seq = seq(1:len);
end

function transformed = apply_forward_transformation(seq, N)
    transformed = [];
    for i = 1:length(seq)
        if seq(i) == 1
            transformed = [transformed, 1, zeros(1, N)];
        else
            transformed = [transformed, 1, zeros(1, N-1)];
        end
    end
end

function transformed = apply_backward_transformation(seq, N)
    transformed = [];
    i = 1;
    while i <= length(seq)
        if seq(i) == 1
            if (i + N <= length(seq) && all(seq(i+1:i+N) == 0))
                transformed = [transformed, 1];
                i = i + N + 1;
            elseif (i + N - 1 <= length(seq) && all(seq(i+1:i+N-1) == 0))
                transformed = [transformed, 0];
                i = i + N;
            else
                return;
            end
        else
            return;
        end
    end
end

function plot_sequence(seq)
    x = 1:25;  % Only the first 25 positions
    y = seq(1:25); % Only the first 25 values
    figure;
    hold on;
    scatter(x(y == 1), ones(sum(y == 1), 1), 50, 'r', 'x', 'LineWidth', 1.5); % Red crosses for 1s
    scatter(x(y == 0), zeros(sum(y == 0), 1), 50, 'b', 'o', 'LineWidth', 1.5); % Blue circles for 0s
    ylim([-0.5 1.5]);
    xlabel('Position in Sequence');
    ylabel('Value');
    hold off;
end
