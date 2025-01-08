#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>

#include <omp.h> // For parallelization (OpenMP)

class optimization {
public:
    /**
     * Simulated Annealing for Optimizing Booster Configurations
     *
     * @param objective_function: A lambda or function that calculates the error given the configuration.
     * @param booster_amounts: A reference to the vector containing initial booster values (modified in-place).
     * @param lower_bound: The lower limit for booster values.
     * @param upper_bound: The upper limit for booster values.
     * @param max_iterations: Maximum number of iterations to perform.
     * @param initial_temperature: Starting temperature for simulated annealing.
     * @param cooling_rate: Rate at which the temperature decreases.
     * @return The best error achieved.
     */
    static double optimize_boosters_simulated_annealing(
        const std::function<double(const std::vector<unsigned short>&)>& objective_function,
        std::vector<unsigned short>& booster_amounts,
        unsigned short lower_bound,
        unsigned short upper_bound,
        int max_iterations = 10000,
        double initial_temperature = 1000.0,
        double cooling_rate = 0.99
    ) {
        // Random number generation setup
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> booster_index_dist(0, booster_amounts.size() - 1);
        std::uniform_int_distribution<int> booster_value_dist(lower_bound, upper_bound);
        std::uniform_real_distribution<double> probability_dist(0.0, 1.0);

        // Initial configuration and evaluation
        std::vector<unsigned short> current_boosters = booster_amounts;
        double current_error = objective_function(current_boosters);

        std::vector<unsigned short> best_boosters = current_boosters;
        double best_error = current_error;

        double temperature = initial_temperature;

        for (int iter = 0; iter < max_iterations; ++iter) {
            // Randomly select a booster to modify
            int index = booster_index_dist(gen);
            unsigned short old_value = current_boosters[index];
            unsigned short new_value = booster_value_dist(gen);

            // Apply the change
            current_boosters[index] = new_value;
            double new_error = objective_function(current_boosters);

            // Accept or reject the new configuration
            if (new_error < current_error || probability_dist(gen) < std::exp((current_error - new_error) / temperature)) {
                current_error = new_error;
                if (new_error < best_error) {
                    best_error = new_error;
                    best_boosters = current_boosters;
                }
            }
            else {
                // Revert the change
                current_boosters[index] = old_value;
            }

            // Decrease temperature
            temperature *= cooling_rate;
        }

        // Update the original booster amounts with the best configuration
        booster_amounts = best_boosters;

        return best_error;
    }

    /**
     * Brute Force Search for Optimizing Booster Configurations
     *
     * @param objective_function: A lambda or function that calculates the error given the configuration.
     * @param booster_amounts: A reference to the vector containing the best booster configuration found (modified in-place).
     * @param lower_bound: The lower limit for booster values.
     * @param upper_bound: The upper limit for booster values.
     * @return The best error achieved.
     */
    static double optimize_boosters_brute_force(
        const std::function<double(const std::vector<unsigned short>&)>& objective_function,
        std::vector<unsigned short>& booster_amounts,
        unsigned short lower_bound,
        unsigned short upper_bound
    ) {
        size_t num_boosters = booster_amounts.size();
        double best_error = std::numeric_limits<double>::infinity();
        std::vector<unsigned short> best_configuration = booster_amounts;

        // Precompute the total number of combinations
        size_t total_combinations = std::pow(upper_bound - lower_bound + 1, num_boosters);
        std::cout << "total_combinations: " << static_cast<long long int>(total_combinations) << std::endl;

        // Parallel brute force using OpenMP
    #pragma omp parallel
        {
            std::vector<unsigned short> current_configuration(num_boosters, lower_bound);
            double local_best_error = std::numeric_limits<double>::infinity();
            std::vector<unsigned short> local_best_configuration = booster_amounts;

    #pragma omp for schedule(dynamic)
            for (long long int combination = 0; combination < static_cast<long long int>(total_combinations); ++combination) {

                // Decode combination into current configuration
                size_t temp = combination;
                for (size_t i = 0; i < num_boosters; ++i) {
                    current_configuration[i] = lower_bound + (temp % (upper_bound - lower_bound + 1));
                    temp /= (upper_bound - lower_bound + 1);
                }

                // Evaluate current configuration
                double current_error = objective_function(current_configuration);
                if (current_error < local_best_error) {
                    local_best_error = current_error;
                    local_best_configuration = current_configuration;

                    // Log the new local configuration
    #pragma omp critical
                    {
                        std::cout << "New local Configuration: ";
                        for (size_t i = 0; i < local_best_configuration.size(); ++i) {
                            std::cout << local_best_configuration[i] << (i < local_best_configuration.size() - 1 ? ", " : " ");
                        }
                        std::cout << " | Error: " << local_best_error << std::endl;
                    }
                }
            }

            // Update global best results and log in a critical section
    #pragma omp critical
            {
                if (local_best_error < best_error) {
                    best_error = local_best_error;
                    best_configuration = local_best_configuration;

                    // Log the new best configuration
                    std::cout << "New Best Configuration: ";
                    for (size_t i = 0; i < best_configuration.size(); ++i) {
                        std::cout << best_configuration[i] << (i < best_configuration.size() - 1 ? ", " : " ");
                    }
                    std::cout << " | Error: " << best_error << std::endl;
                }
            }
        }

        // Update the booster amounts with the best configuration
        booster_amounts = best_configuration;

        return best_error;
    }


    /**
     * Particle Swarm Optimization for Optimizing Booster Configurations
     *
     * @param objective_function: A lambda or function that calculates the error given the configuration.
     * @param booster_amounts: A reference to the vector containing the best booster configuration found (modified in-place).
     * @param lower_bound: The lower limit for booster values.
     * @param upper_bound: The upper limit for booster values.
     * @param swarm_size: Number of particles in the swarm.
     * @param max_iterations: Maximum number of iterations to perform.
     * @param inertia_weight: Weight for the particle's previous velocity.
     * @param cognitive_weight: Weight for the particle's best-known position.
     * @param social_weight: Weight for the swarm's global best position.
     * @return The best error achieved.
     */
    static double optimize_boosters_particle_swarm(
        const std::function<double(const std::vector<unsigned short>&)>& objective_function,
        std::vector<unsigned short>& booster_amounts,
        unsigned short lower_bound,
        unsigned short upper_bound,
        int swarm_size = 30,
        int max_iterations = 1000,
        double inertia_weight = 0.5,
        double cognitive_weight = 1.5,
        double social_weight = 1.5
    ) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist_real(0.0, 1.0);
        std::uniform_int_distribution<int> dist_int(lower_bound, upper_bound);

        size_t num_boosters = booster_amounts.size();

        // Initialize particles
        std::vector<std::vector<unsigned short>> particles(swarm_size, std::vector<unsigned short>(num_boosters));
        std::vector<std::vector<unsigned short>> particle_best_positions = particles;
        std::vector<double> particle_best_scores(swarm_size, std::numeric_limits<double>::infinity());

        std::vector<unsigned short> global_best_position = booster_amounts;
        double global_best_score = std::numeric_limits<double>::infinity();

        std::vector<std::vector<double>> velocities(swarm_size, std::vector<double>(num_boosters, 0.0));

        // Initialize particle positions and velocities
        for (int i = 0; i < swarm_size; ++i) {
            for (size_t j = 0; j < num_boosters; ++j) {
                particles[i][j] = dist_int(gen);
                velocities[i][j] = dist_real(gen) * (upper_bound - lower_bound) / 2.0;
            }
            particle_best_positions[i] = particles[i];
            particle_best_scores[i] = objective_function(particles[i]);

            if (particle_best_scores[i] < global_best_score) {
                global_best_score = particle_best_scores[i];
                global_best_position = particle_best_positions[i];
            }
        }

        // PSO main loop
        for (int iter = 0; iter < max_iterations; ++iter) {
            for (int i = 0; i < swarm_size; ++i) {
                for (size_t j = 0; j < num_boosters; ++j) {
                    // Update velocity
                    velocities[i][j] = inertia_weight * velocities[i][j]
                        + cognitive_weight * dist_real(gen) * (particle_best_positions[i][j] - particles[i][j])
                        + social_weight * dist_real(gen) * (global_best_position[j] - particles[i][j]);

                    // Clamp velocity and update position
                    velocities[i][j] = std::max(std::min(velocities[i][j], (double)upper_bound), (double)lower_bound);
                    particles[i][j] = std::max(std::min((int)(particles[i][j] + velocities[i][j]), (int)upper_bound), (int)lower_bound);
                }

                // Evaluate particle
                double score = objective_function(particles[i]);
                if (score < particle_best_scores[i]) {
                    particle_best_scores[i] = score;
                    particle_best_positions[i] = particles[i];
                }

                // Update global best
                if (score < global_best_score) {
                    global_best_score = score;
                    global_best_position = particles[i];
                }
            }
        }

        // Update the booster amounts with the best configuration
        booster_amounts = global_best_position;

        return global_best_score;
    }

    /**
     * Genetic Algorithm for Optimizing Booster Configurations
     *
     * @param objective_function: A lambda or function that calculates the error given the configuration.
     * @param booster_amounts: A reference to the vector containing the best booster configuration found (modified in-place).
     * @param lower_bound: The lower limit for booster values.
     * @param upper_bound: The upper limit for booster values.
     * @param population_size: Number of individuals in the population.
     * @param generations: Maximum number of generations to evolve.
     * @param mutation_rate: Probability of mutating a single booster value.
     * @param crossover_rate: Probability of performing crossover between two parents.
     * @return The best error achieved.
     */
    static double optimize_boosters_genetic(
        const std::function<double(const std::vector<unsigned short>&)>& objective_function,
        std::vector<unsigned short>& booster_amounts,
        unsigned short lower_bound,
        unsigned short upper_bound,
        int population_size = 50,
        int generations = 100,
        double mutation_rate = 0.1,
        double crossover_rate = 0.7
    ) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist_real(0.0, 1.0);
        std::uniform_int_distribution<int> dist_int(lower_bound, upper_bound);

        size_t num_boosters = booster_amounts.size();

        // Initialize population
        std::vector<std::vector<unsigned short>> population(population_size, std::vector<unsigned short>(num_boosters));
        for (auto& individual : population) {
            for (auto& booster : individual) {
                booster = dist_int(gen);
            }
        }

        // Evaluate initial population
        std::vector<double> fitness(population_size);
        for (int i = 0; i < population_size; ++i) {
            fitness[i] = objective_function(population[i]);
        }

        // Main genetic algorithm loop
        for (int generation = 0; generation < generations; ++generation) {
            // Selection: Sort population by fitness
            std::vector<size_t> sorted_indices(population_size);
            std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
            std::sort(sorted_indices.begin(), sorted_indices.end(),
                [&fitness](size_t a, size_t b) { return fitness[a] < fitness[b]; });

            // Elitism: Keep the best individual
            std::vector<std::vector<unsigned short>> new_population;
            new_population.push_back(population[sorted_indices[0]]);

            // Generate offspring
            while (new_population.size() < population_size) {
                // Select parents
                size_t parent1_idx = sorted_indices[rand() % (population_size / 2)];
                size_t parent2_idx = sorted_indices[rand() % (population_size / 2)];

                // Crossover
                std::vector<unsigned short> offspring = population[parent1_idx];
                if (dist_real(gen) < crossover_rate) {
                    for (size_t i = 0; i < num_boosters; ++i) {
                        if (dist_real(gen) < 0.5) {
                            offspring[i] = population[parent2_idx][i];
                        }
                    }
                }
                
                // Mutation
                for (size_t i = 0; i < num_boosters; ++i) {
                    if (dist_real(gen) < mutation_rate) {
                        offspring[i] = dist_int(gen);
                    }
                }

                new_population.push_back(offspring);
            }

            // Replace old population with new one
            population = new_population;

            // Recompute fitness
            for (int i = 0; i < population_size; ++i) {
                fitness[i] = objective_function(population[i]);
            }
        }

        // Return the best individual0
        auto best_idx = std::min_element(fitness.begin(), fitness.end()) - fitness.begin();
        booster_amounts = population[best_idx];
        return fitness[best_idx];
    }

};

#endif // OPTIMIZATION_H
