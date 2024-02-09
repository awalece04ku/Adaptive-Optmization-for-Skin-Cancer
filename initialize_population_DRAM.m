function population = initialize_population_DRAM(n, d, lb, ub,fitness)
%% implementation of the Delayed-Rejection Adaptive Metropolis (DRAM) method in MATLAB:
    % n - number of individuals in the population
    % d - number of dimensions (variables) of each individual
    % lb - lower bound of each dimension
    % ub - upper bound of each dimension

    population = zeros(n, d);

    % Initialize first individual randomly
    population(1, :) = lb + (ub - lb) .* rand(1, d);

    for i = 2:n
        % Generate candidate individual
        candidate = lb + (ub - lb) .* rand(1, d);

        % Calculate acceptance probability
        acceptance_prob = min(1, fitness(candidate) / fitness(population(i-1, :)));

        % Accept or reject candidate with probability
        if rand() <= acceptance_prob
            population(i, :) = candidate;
        else
            % Generate second candidate
            candidate2 = lb + (ub - lb) .* rand(1, d);

            % Calculate second acceptance probability
            acceptance_prob2 = min(1, fitness(candidate2) / fitness(population(i-1, :)));

            % Accept or reject second candidate with probability
            if rand() <= acceptance_prob2
                population(i, :) = candidate2;
            else
                population(i, :) = population(i-1, :);
            end
        end
    end
end
