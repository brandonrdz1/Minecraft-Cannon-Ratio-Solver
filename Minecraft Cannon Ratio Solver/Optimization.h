#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>

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
        double cooling_rate = 0.999999
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
};

#endif // OPTIMIZATION_H
