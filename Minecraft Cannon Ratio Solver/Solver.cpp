#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>

#include "optimization.h"

#define center -132111.5 // Where everything starts from
#define center_x -132111.5 // Where everything starts from

namespace hammer {
	double explosion_height = 0.061250001192092896;
	double gravity = -0.04;
	double drag = 0.98;

	void compute_booster_to_power_x_exposure(double power_x0, double* booster_x_values_G, double* booster_to_power_x_exposures) {
		std::cout << "Booster to Hammer's Power Exposure: " << std::endl;
		for (int i = 0; i < 12; i++) {
			double distance = power_x0 - booster_x_values_G[i];
			double f = (1.0 - distance / 8.0) * 1.0 / distance; // Just to keep correct calculations TODO: simplify
			booster_to_power_x_exposures[i] = distance * f;			// Output

			std::cout << "  booster[" << i << "] (blocks/tick): " << booster_to_power_x_exposures[i] << std::endl;
		}
		return;
	}

	double position_after_17gt_L(double* booster_to_power_x_exposures, unsigned short int* booster_amounts, double power_x0_L, unsigned int power_amount) {
		// X Position of Hammer's Power
		double power_x1_L = power_x0_L;
		// Affected by boosters
		for (int i = 0; i < 12; i++) {
			power_x1_L += booster_amounts[i] * booster_to_power_x_exposures[i];
		}

		double hammer_x0_L = 0.0;							// Local initial position of the hammer
		double distance = hammer_x0_L - power_x1_L;		// Distance between Hammer's power and Hammer
		double f = (1.0 - distance / 8.0) * 1.0 / distance; // Force

		// Freefall of the hammer for 17 gt
		double hammer_u0 = f * power_amount * distance;	// Initial velocity of the Hammer
		double hammer_x17_L = hammer_x0_L;
		for (int i = 0; i < 17; i++) { // 17
			hammer_x17_L += hammer_u0;
			hammer_u0 *= 0.98;
		}

		// Return's Hammer's local x position
		return hammer_x17_L;
	}

	void compute_booster_x_Hrange_virtual_exposure(unsigned short int* booster_amounts, double* booster_toHP_x_exposures, double power_H_x1i_L, unsigned int power_H_amount, double* booster_x_Hrange_virtual_exposure) {
		// Compute the reference position
		std::cout << "Booster to Hammer's Final Position Exposure: " << std::endl;
		for (int i = 0; i < 12; i++) { booster_amounts[i] = 6; }
		double hammer_x_ref_final_L = position_after_17gt_L(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, power_H_amount);

		// Find the effects of pertubations
		for (int i = 0; i < 12; i++) {
		// Decrement 1
		booster_amounts[i] = 5;
		double hammer_x_neg_perturbation_final_L = position_after_17gt_L(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, power_H_amount);
		double decrement_rel_change = hammer_x_ref_final_L - hammer_x_neg_perturbation_final_L;

		// Increment 1
		booster_amounts[i] = 7;
		double hammer_x_pos_perturbation_final_L = position_after_17gt_L(booster_toHP_x_exposures, booster_amounts, power_H_x1i_L, power_H_amount);
		double increment_rel_change = hammer_x_pos_perturbation_final_L - hammer_x_ref_final_L;

		// Compute Inderect relation of incrementing / decrementing the booster[i]
		booster_x_Hrange_virtual_exposure[i] = (decrement_rel_change + increment_rel_change) / 2.0;
		// Print out virtual exposure
		std::cout << "  booster[" << i << "] (blocks): " << booster_x_Hrange_virtual_exposure[i] << std::endl;

		// reset
		booster_amounts[i] = 6;
		}

	}

	double compute_virtual_initial_position(double* booster_x_range_virtual_exposure, double x_ref_final_L) {
		double ref_booster_dist = 0.0;
		for (int i = 0; i < 12; i++) {
			ref_booster_dist += 6 * booster_x_range_virtual_exposure[i];
		}
		double x0virtual_L = x_ref_final_L - ref_booster_dist;
		std::cout << "Initial Virtual Hammer's Position (local): " << x0virtual_L << std::endl;
		return x0virtual_L;
	}

	double compute_objective(
		double* booster_toHP_x_exposures,			// To check if the Hammer's power is in range
		unsigned short int* booster_amounts,		// To check if the Hammer's power is in range
		double power_H_x1i_L,						// To check if the Hammer's power is in range
		double* booster_x_Hrange_virtual_exposure,	// To compute Hammer's final position
		double hammer_x0virtual_L,					// Hammer's virtual initial position
		double target_x_L							// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_H_x1_L = power_H_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_H_x1_L += booster_amounts[i] * booster_toHP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_H_x1_L > -0.5 || power_H_x1_L < -7.5) { // TODO: Play with bounds
			return 10000000.0; // should be infinity
		}
		// Find the final position of the Hammer
		double hammer_x17_L = hammer_x0virtual_L;
		for (int i = 0; i < 12; i++) {
			hammer_x17_L += booster_amounts[i] * booster_x_Hrange_virtual_exposure[i];
		}
		// Compute the Error
		double error = abs(hammer_x17_L - target_x_L);
		return error;
	}

	double debug_objective(
		double* booster_toHP_x_exposures,				// To check if the Hammer's power is in range
		std::vector<unsigned short>& booster_amounts,	// To check if the Hammer's power is in range
		double power_H_x1i_L,							// To check if the Hammer's power is in range
		double* booster_x_Hrange_virtual_exposure,		// To compute Hammer's final position
		double hammer_x0virtual_L,						// Hammer's virtual initial position
		double target_x_L								// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_H_x1_L = power_H_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_H_x1_L += booster_amounts[i] * booster_toHP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_H_x1_L > -0.5 || power_H_x1_L < -7.5) { // TODO: Play with bounds
			// return invalid
			// std::cout << "  Hammer's power is out of range/ahead of it" << std::endl; // Penalty funciton
			return 10000000.0; // should be infinity
		}
		// else
		std::cout << "  Hammer's Power x position (local): " << power_H_x1_L << std::endl;
		std::cout << "  Hammer's Power x position (global): " << power_H_x1_L + center << std::endl;
		// Find the final position of the Hammer
		double hammer_x17_L = hammer_x0virtual_L;
		for (int i = 0; i < 12; i++) {
			hammer_x17_L += booster_amounts[i] * booster_x_Hrange_virtual_exposure[i];
		}
		std::cout << "  Hammer at 17gt x position (local): " << hammer_x17_L << std::endl;
		std::cout << "  Hammer at 17gt x position (global): " << hammer_x17_L + center << std::endl;
		// Compute the Error
		double error = abs(hammer_x17_L - target_x_L);
		return error;
	}
}

namespace arrow {

	double arrow_eye = 0.12999999523162842;
	double gravity = -0.05f;
	double drag = 0.99f;

	void compute_booster_to_power_x_exposure(double power_x0, double* booster_x_values_G, double* booster_to_power_x_exposures) {
		std::cout << "Booster to Arrows's Power Exposure: " << std::endl;
		for (int i = 0; i < 12; i++) {
			double distance = power_x0 - booster_x_values_G[i];
			double f = (1.0 - distance / 8.0) * 1.0 / distance; // Just to keep correct calculations TODO: simplify
			booster_to_power_x_exposures[i] = distance * f;			// Output

			std::cout << "  booster[" << i << "] (blocks/tick): " << booster_to_power_x_exposures[i] << std::endl;
		}
		return;
	}

	void compute_arrow_x17gt_L_map(double* arrow_x17_map, double* arrow_u17_map, double power_A_amount) {
		// double x_min = -7.5, x_max = -0.5, step = 0.05;
		// double power_x_arr[141];
		// for (int i = 0; i <= 140; ++i) { power_x_arr[i] = x_min + i * step; }
		// for (int i = 0; i <= 140; ++i) { std::cout << "power_loc_arr[" << i << "] = " << power_x_arr[i] << std::endl; }
		
		double dy_squared = arrow_eye * arrow_eye;
		for (int i = 0; i <= 140; i++) {
			double dx = -7.5 + i * 0.05;
			double distance = std::sqrt(dx * dx + dy_squared);
			double f = (1.0 - distance / 8.0);
			double f0 = f;

			// elaborate way to compute initial velocity, similar to game
			for (int j = 1; j < power_A_amount; j++) { f0 += f; }
			
			// freefall
			double x = 0.0;		// initial local position of the arrow
			double v = f0;		// initial velocity of the arrow
			for (int j = 0; j < 17; j++) {
				x += v;
				v *= drag;
			}
			arrow_x17_map[i] = x;
			arrow_u17_map[i] = v;

			//std::cout << dx << ", " << arrow_u17_map[i] << std::endl;
		}
		return;
	}

	double compute_objective(
		double* booster_toAP_x_exposures,			// To check if the Hammer's power is in range
		unsigned short int* booster_amounts,		// To check if the Hammer's power is in range
		double power_A_x1i_L,						// To check if the Hammer's power is in range
		double* arrow_x17_map,	
		double target_x_L							// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_A_x1_L = power_A_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_A_x1_L += booster_amounts[i] * booster_toAP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_A_x1_L > -1 || power_A_x1_L < -7.5) { // TODO: Play with bounds
			return 10000000.0; // should be infinity
		}
		
		// Map Lookup
		// Find the final position of the Hammer
		int lower_index = static_cast<int>((power_A_x1_L - -7.5) / 0.05);
		int upper_index = lower_index + 1;

		// Fetch corresponding values from the map
		double lower_value = arrow_x17_map[lower_index];
		double upper_value = arrow_x17_map[upper_index];

		// Linear interpolation
		double fraction = (power_A_x1_L - (-7.5 + lower_index * 0.05)) / 0.05;
		double arrow_x17_L = lower_value + fraction * (upper_value - lower_value);

		// Compute the error
		double error = std::abs(arrow_x17_L - target_x_L);
		return error;
	}

	double debug_objective(
		double* booster_toAP_x_exposures,				// To check if the Hammer's power is in range
		std::vector<unsigned short>& booster_amounts,	// To check if the Hammer's power is in range
		double power_A_x1i_L,							// To check if the Hammer's power is in range
		double* arrow_x17_map,
		double target_x_L								// To minimize
	) {
		// Check if this will place the hammer's power within range
		// Position of hammer's power with current boosters:
		double power_A_x1_L = power_A_x1i_L;
		for (int i = 0; i < 12; i++) {
			power_A_x1_L += booster_amounts[i] * booster_toAP_x_exposures[i];
		}
		// Check if the Hammer's power is within explosion range
		if (power_A_x1_L > -1 || power_A_x1_L < -7.5) { // TODO: Play with bounds
			return 10000000.0; // should be infinity
		}

		// Arrow's Power initial position
		std::cout << "  Arrow's Power x position (local): " << power_A_x1_L << std::endl;
		std::cout << "  Arrow's Power x position (global): " << power_A_x1_L + center << std::endl;

		// Map Lookup
		// Find the final position of the Hammer
		int lower_index = static_cast<int>((power_A_x1_L - -7.5) / 0.05);
		int upper_index = lower_index + 1;

		// Fetch corresponding values from the map
		double lower_value = arrow_x17_map[lower_index];
		double upper_value = arrow_x17_map[upper_index];

		// Linear interpolation
		double fraction = (power_A_x1_L - (-7.5 + lower_index * 0.05)) / 0.05;
		double arrow_x17_L = lower_value + fraction * (upper_value - lower_value);
		std::cout << "  Arrow at 17gt x position (local): " << arrow_x17_L << std::endl;
		std::cout << "  Arrow at 17gt x position (global): " << arrow_x17_L + center << std::endl;

		// Compute the error
		double error = std::abs(arrow_x17_L - target_x_L);
		return error;
	}
}

class Tnt {
public:
	double x = 0.0, y = 0.0, z = 0.0;
	double u = 0.0, v = 0.0, w = 0.0;
	unsigned int amount = 1;

	static const double explosion_height, gravity, drag;
	static const double distance(const Tnt& a, const Tnt& b) {
		double dx = a.x - b.x;
		double dy = a.y - b.y;
		double dz = a.z - b.z;
		return std::sqrt(dx * dx + dy * dy + dz * dz);
	}

	Tnt(const std::string& attributes, unsigned int explosivePower = 1.0) : amount(explosivePower) {
		std::istringstream stream(attributes);
		stream >> x >> y >> z >> u >> v >> w;
	};

	Tnt(double x, double y, double z, double u = 0.0, double v = 0.0, double w = 0.0, unsigned int explosivePower = 1.0) : x(x), y(y), z(z), u(u), v(v), w(w), amount(explosivePower) {};

	void explosion(Tnt source, double exposure = 1.0) {
		double dx = x - source.x;
		double dy = y - source.y; // double dy = y - (source.y + explosion_height);
		double dz = z - source.z;
		float distance = std::sqrt(dx * dx + dy * dy + dz * dz); // double
		double f = ( -1.0 / 8.0 + 1.0 / distance) * source.amount * exposure;
		u += f * dx;
		v += f * dy;
		w += f * dz;
	}

	void swing(Tnt source) {
		Tnt origin(source);
		std::vector<std::vector<double>> locations;
		double dx, dy, dz;
		float distance;
		double f;
		for (int i = 0; i < source.amount; i++) {
			Tnt next = origin;
			next.freefall(1);
			if (next.y > 254.0199999809265) { next.y = 254.0199999809265; }
			std::vector<double> loc(3);
			loc[0] = next.x;
			loc[1] = next.y;
			loc[2] = next.z;
			locations.push_back(loc);
			dx = origin.x - locations[i][0], dy = origin.y - locations[i][1], dz = origin.z - locations[i][2];
			distance = std::sqrt(dx * dx + dy * dy + dz * dz);
			if (distance > 8.0 || distance == 0.0) {
				continue;
			}
			f = -1.0/8.0 + 1.0 / distance;
			origin.u += dx * f;
			origin.v += dy * f;
			origin.w += dz * f;
		}

		for (int i = 0; i < locations.size(); i++) {
			Tnt source(locations[i][0], locations[i][1], locations[i][2]);
			//source.print("Source[" + std::to_string(i) + "]: ");
			this->explosion(source);
		}
	}

	void freefall(int ticks) {
		// print(" 0 ");
		for (int i = 0; i < ticks; i++) {
			v += Tnt::gravity;

			x += u;
			y += v;
			z += w;
			// print(" " + std::to_string(i + 1) + " ");
			u *= Tnt::drag;
			v *= Tnt::drag;
			w *= Tnt::drag;
		}
	}

	std::vector<double> getVel() {
		return { u, v, w };
	}

	void print(std::string str = "") {
		std::cout << str << "x: " << x << "\ty: " << y << "\tz: " << z << "\t| u: " << u << "\tv: " << v << "\tw: " << w << std::endl;
	}
};
// Define static const members outside the class
const double Tnt::explosion_height = 0.061250001192092896;
const double Tnt::gravity = -0.04f;
const double Tnt::drag = 0.98f;

double CalculateTargetError(double ypos) {
	double targetError = std::numeric_limits<double>::infinity(); // Best error found so far
	double low = 0.0;
	double high = 10.1;
	double targetY = 255.0;

	double bestVelocity = low; // Best velocity corresponding to the minimum error
	std::cout << std::endl;
	// Create binary search with robust error handling
	while (high - low > 1e-16) { // Use a more practical precision threshold
		double mid = (low + high) / 2.0;
		Tnt testTnt(0, ypos, 0, 0, 0, 0); // Reset Tnt object for each test
		testTnt.v = mid;
		testTnt.freefall(7);

		double error = std::abs(testTnt.y - targetY); // Calculate the current error

		// Update the best error and velocity if the current one is better
		if (error < targetError) {
			targetError = error;
			bestVelocity = mid;
		}

		// Adjust bounds based on the position of the current error
		if (testTnt.y < targetY) {
			low = mid; // Test Tnt falls short, increase velocity
		}
		else {
			high = mid; // Test Tnt overshoots, decrease velocity
		}

		// Debugging output
		std::cout << "mid: " << mid
			<< "\terror: " << error
			<< "\ty: " << testTnt.y
			<< "\ttargetError: " << targetError << std::endl;

		// Add a fallback termination condition
		if (targetError < 1e-16) { // Acceptable error threshold
			break;
		}
	}

	// Final output
	std::cout << "Best Velocity: " << bestVelocity
		<< " with Error: " << targetError << std::endl;
	std::cout << std::endl;
	return bestVelocity;
}

bool generateNextPermutation(int currentArray[], const int startArray[], const int finalArray[], size_t size) {
	size_t i = size;

	while (i > 0) {
		--i;
		if (currentArray[i] < finalArray[i]) {
			++currentArray[i];
			for (size_t j = i + 1; j < size; ++j) {
				currentArray[j] = startArray[j]; // Reset subsequent values to their start
			}
			return true;
		}
	}
	return false;
}

int main() {
	/************* Output Formatting *************/
	std::cout << std::setprecision(17);

	/* In-game Conditions (Global) */
	//// Positions/Velocities
	Tnt swingTnt("-199971.5000002721 254.0199999809265 -333113.50999999046 0.9701997615808802 0.0 0.0", 75);
	Tnt swingTntRev("-199971.5 250.0 -333113.5 0.0 0.0 0.0", 1);
	swingTnt.explosion(swingTntRev); // keeps it up
	Tnt powerTnt0("-199971.49755261914 254.0199999809265 -333108.5534764172 -0.00740155756734615 0.0 0.9177931383509221");
	Tnt powerTnt0Rev0("-199964.5 254.0 -333111.5", 1);
	Tnt powerTnt0Rev1("-199978.5 254.0 -333111.5", 1);
	powerTnt0.explosion(powerTnt0Rev0, (float)1.0/3.0);
	powerTnt0.explosion(powerTnt0Rev1, (float)1.0/3.0);
	powerTnt0.print("Before Swing: ");

	int count = 0;
	int add = -1;
	std::vector<double> xlocations;
	std::vector<double> ylocations;
	for (int i = 75; i <= 165; i++) {
		add++;
		Tnt powerTnt1 = powerTnt0; // reset
		swingTnt.amount = i;
		powerTnt1.swing(swingTnt);
		//powerTnt1.print("After Swing NoFF with [" + std::to_string(i) + "] power: ");
		
		powerTnt1.freefall(1);
		double decimal = powerTnt1.x - std::floor(powerTnt1.x);
		
		if (decimal < 0.49f || decimal > 0.51f) {
			
			continue;
		}
		count++;
		// guider alignment
		// powerTnt1.y = 254.0; TODO: check if this is necessary
		powerTnt1.z = -333100.49000000954; // guider coordinate

		powerTnt1.print("After Swing + FF with [" + std::to_string(i) + "] (+" +std::to_string(add) +") power: ");
		xlocations.push_back(powerTnt1.x);
		ylocations.push_back(powerTnt1.y);
		add = 0;
	}
	// we use 0-14 but index = 4 has a maximum of 0 power

	std::cout << "count: " << count << std::endl;
	// generate micro boosters
	Tnt* microBoosters[15];
	for (int i = 0; i < 15; i++) {
		microBoosters[i] = new Tnt(xlocations[i], ylocations[i], -333382.49000000954, 0.0, 0.0, 0.0, 1); // create the micro boosters array with 1 tnt each
		microBoosters[i]->print("MicroBooster[" + std::to_string(i) + "]:\t");
	}
	// start the projectile
	Tnt projectile0("-199971.50999999046 253.97999998182058 -333382.49000000954 -1.9164525895466311 -0.03919999988675116 0.0");
	//Tnt projectile0("-199973.46103707515 253.97999998182058 -333382.49000000954 -1.9164525895466311 -0.03919999988675116 0.0");
	// compute the boosters to power exposure
	std::vector<std::vector<double>> microExposures(15);
	for (int i = 0; i < 15; i++) {
		Tnt projectile1i(projectile0.x, projectile0.y, projectile0.z, 0.0, 0.0, 0.0); // dummy inertial projectile to get exposure effects
		projectile1i.explosion(*microBoosters[i]);
		microExposures[i] = projectile1i.getVel();
		std::cout << "microBooster["<<i<<"] effects: " << microExposures[i][0] << ", " << microExposures[i][1] << ", " << microExposures[i][2] << std::endl;
		double dist = Tnt::distance(projectile0, *microBoosters[i]);
	}

	///// debugging
	double min_power = 0;
	double max_power = 7.53333;
	double mid_power = max_power / 2.0;

	// compute min/max influence
	double minX_influence = 0.0;
	double maxX_influence = 0.0;
	double midX_influence = 0.0;
	double minY_influence = 0.0;
	double maxY_influence = 0.0;
	double midY_influence = 0.0;

	for (int i = 0; i < 15; i++) {
		minX_influence += min_power * microExposures[i][0];
		minY_influence += min_power * microExposures[i][1];
		midX_influence += mid_power * microExposures[i][0];
		
		maxX_influence += max_power * microExposures[i][0];
		maxY_influence += max_power * microExposures[i][1];
		midY_influence += mid_power * microExposures[i][1];
	}

	std::cout << "minX_influence: " << minX_influence << std::endl;
	std::cout << "maxX_influence: " << maxX_influence << std::endl;
	std::cout << "midX_influence: " << midX_influence << std::endl;

	std::cout << "minY_influence: " << minY_influence << std::endl;
	std::cout << "maxY_influence: " << maxY_influence << std::endl;
	std::cout << "midY_influence: " << midY_influence << std::endl;
	
	// find the taget yvelocity
	// Find the target y-velocity
	
	double targetYvel0 = CalculateTargetError(projectile0.y);
	std::cout << "targetYvel0: " << targetYvel0 << std::endl;

	// add the influence of the ypower
	Tnt yPower0("-199971.49000000954 253.97730819807788 -333382.49000000954 0.0 0.0 0.0");
	yPower0.print("yPower0: ");
	
	int record = 0;
	double errorYvel = 100000000.0;
	double startVel255;
	for (int i = 0; i < 10000; i++) {
		yPower0.amount = i;
		Tnt projectile1 = projectile0; // reset

		// projectile1.u += midX_influence;
		// projectile1.v += midY_influence;

		projectile1.explosion(yPower0);
		double errorCalculated = std::abs(projectile1.v - targetYvel0);
		if (errorCalculated < errorYvel && projectile1.v > targetYvel0) {
			errorYvel = errorCalculated;
			record = i;
			std::cout << "record: " << record << std::endl;
			std::cout << "errorYvel: " << errorYvel << std::endl;

			std::cout << "xvel: " << projectile1.u << std::endl;
			std::cout << "yvel: " << projectile1.v << std::endl;
			startVel255 = projectile1.v;
			// history:

			Tnt projectile2 = projectile1;
			projectile2.print("Projectile Freefall: ");
			for (int j = 0; j < 15; j++) {
				projectile2.freefall(1);
				projectile2.print();
			}
			std::cout << std::endl;
		}
	}
	//////////////////
	/*
	Tnt proj0("-199971.49000000954 226.0 -333479.49000000954 0.0 0.0 0.0");
	Tnt bo("-199971.49000000954 225.0 -333480.51 0.0 0.0 0.0");
	Tnt ss0("-199971.49000000954 225.875 -333482.5");
	Tnt ss1("-199971.49000000954 225.875 -333484.5099");

	for (int i = 0; i < 100; i++) {
		for (int j = 0; j < 20; j++) {
			for (int k = 0; k < 1; k++) {

				bo.amount = i;
				ss0.amount = j;
				ss1.amount = k;

				Tnt proj1 = proj0;

				proj1.explosion(bo);
				proj1.explosion(ss0);
				proj1.explosion(ss1);

				proj1.freefall(1);

				if (proj1.y > 253.9774 - 0.002 && proj1.y < 253.9774 + 0.0001) {
					std::cout << "Amount: " << i << " " << j << " " << k << std::endl;
					std::cout << " Height: " << proj1.y << std::endl;
				}

			}
			
		}
	}*/

	// 255.0 script
	int microAmountStart[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int microAmountCurrent[12] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	int microAmountFinal[12] = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 };

	double finalVelocity255 = 0.0;

	// startVel255 contains the starting velocity of the projectile
	std::cout << "starting velocity: " << startVel255 << std::endl;
	std::cout << "starting error: " << startVel255 - targetYvel0 << std::endl;

	Tnt rangeOutput(0.0, projectile0.y, 0.0, 0.0, startVel255, 0.0);
	rangeOutput.freefall(7);
	std::cout << "minimum possible value: " << rangeOutput.y << std::endl;
	for (int i = 0; i < 12; i++) {
		finalVelocity255 += (12 * microExposures[i][1]);
	}
	rangeOutput = Tnt(0.0, projectile0.y, 0.0, 0.0, startVel255 + finalVelocity255, 0.0);
	rangeOutput.freefall(7);
	std::cout << "maximum possible value: " << rangeOutput.y << std::endl;


	double record255;
	double error255 = 100000000.0;
	do {
		finalVelocity255 = startVel255;
		for (int i = 0; i < 12; ++i) {
			finalVelocity255 += (microAmountCurrent[i] * microExposures[i][1]);
		}
		if (std::abs(finalVelocity255 - targetYvel0) <= error255) {
			for (int i = 0; i < 12; i++) {
				std::cout << microAmountCurrent[i] << " ";
			} 
			std::cout << std::endl;

			error255 = std::abs(finalVelocity255 - targetYvel0);
			record255 = finalVelocity255;
			std::cout << "  record255: " << record255 << std::endl;
			std::cout << "  error255: " << error255 << std::endl;

			Tnt output(0.0, projectile0.y, 0.0, 0.0, finalVelocity255, 0.0);
			output.freefall(7);
			std::cout << "  final position: " << output.y << std::endl;
		}
	} while (generateNextPermutation(microAmountCurrent, microAmountStart, microAmountFinal, 12));

	for (int i = 0; i < 15; i++) {
		delete microBoosters[i];
	}
	return 0;
}

