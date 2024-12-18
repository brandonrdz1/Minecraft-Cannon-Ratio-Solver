#include <iostream>
#include <iomanip>
#include <cmath>

double objective_function(double* booster_exposures_x, unsigned short int* booster_amounts, double power_H_x0, unsigned int power_H_amount, double target_x) {
	// not yet goning to return minimum or max
	double power_H_x1 = power_H_x0;
	for (int i = 0; i < 12; i++) {
		power_H_x1 += booster_amounts[i] * booster_exposures_x[i];
	}
	double distance = 0.0 - power_H_x1; // power_H to Hammer distance
	double f = (1.0 - distance / 8.0) * 1.0 / distance;
	
	double hammer_u0 = f * power_H_amount * distance;
	double hammer_x1 = hammer_u0;
	// std::cout << "First position of hammer x Local: " << hammer_x1 << std::endl;

	return hammer_x1;
}

int main() {

	std::cout << std::setprecision(17);
	//// In-game Conditions
	double center_x = -86214.5;						//  -86214.587930.5

	double booster_x_values[12] = { 
		-86234.78586161736,
		-86235.03129429705,
		-86235.29123018472,
		-86235.56509943453,
		-86235.85226417152,
		-86236.15201904518,
		-86233.91205295271,
		-86234.25911295973,
		-86234.62067852802,
		-86234.99617983929,
		-86235.38497905637,
		-86235.78637088025
	};

	// State before afected by ratio and constant booster
	double power_H_x0_value = -86228.49000000954;	// -86228.49000000954 278.06125000119226 87930.5
	double power_H_u0_value = 0.0;					// 0.0 13.581624504779418 -2.4523938435550006E-15

	// State after constant booster only
	double power_H_x1_value = -86232.08075066289;	// -86232.08075066289 316.3339737459176 87930.5
	double power_H_u1_value = -3.5189356402947953;		// -3.5189356402947953 37.507269269830815 -1.158739770801276E-14
	
	// Function parameters
	double booster_x_exposures[12]; // Is also the distance it pushes the hammer along x
	std::cout << "Booster Exposure: " << std::endl;
	for (int i = 0; i < 12; i++) {
		double distance = power_H_x0_value - booster_x_values[i];
		double f = (1.0 - distance / 8.0) * 1.0 / distance;
		booster_x_exposures[i] = distance * f;
		std::cout << "  booster[" << i << "]: " << booster_x_exposures[i] << std::endl;
	}
	double power_H_x0 = power_H_x1_value - center_x; // Virtual Starting location of the hammer
	unsigned int power_H_amount = 36;

	// Function Target
	double target_x = 100.0;


	
	// Result?
	double hammer_x1;
	unsigned short int booster_amount[12] = { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }; // Limited to Between 0-12
	
	std::cout << "{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 } Booster Configuration: " << std::endl;
	hammer_x1 = objective_function(booster_x_exposures, booster_amount, power_H_x0, power_H_amount, target_x) + center_x;
	std::cout << "  First position of hammer x Global: " << hammer_x1 << std::endl;

	booster_amount[11] = 7;
	std::cout << "{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7 } Booster Configuration: " << std::endl;
	hammer_x1 = objective_function(booster_x_exposures, booster_amount, power_H_x0, power_H_amount, target_x) + center_x;
	std::cout << "  First position of hammer x Global: " << hammer_x1 << std::endl;

	booster_amount[11] = 6;
	booster_amount[10] = 7;
	std::cout << "{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 6 } Booster Configuration: " << std::endl;
	hammer_x1 = objective_function(booster_x_exposures, booster_amount, power_H_x0, power_H_amount, target_x) + center_x;
	std::cout << "  First position of hammer x Global: " << hammer_x1 << std::endl;

	booster_amount[11] = 7;
	std::cout << "{ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7 } Booster Configuration: " << std::endl;
	hammer_x1 = objective_function(booster_x_exposures, booster_amount, power_H_x0, power_H_amount, target_x) + center_x;
	std::cout << "  First position of hammer x Global: " << hammer_x1 << std::endl;

	booster_amount[0] = 0;
	booster_amount[1] = 0;
	booster_amount[11] = 11;
	std::cout << "{ 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 7, 11 } Booster Configuration: " << std::endl;
	hammer_x1 = objective_function(booster_x_exposures, booster_amount, power_H_x0, power_H_amount, target_x) + center_x;
	std::cout << "  First position of hammer x Global: " << hammer_x1 << std::endl;

	booster_amount[11] = 12;
	std::cout << "{ 0, 0, 6, 6, 6, 6, 6, 6, 6, 6, 7, 12 } Booster Configuration: " << std::endl;
	hammer_x1 = objective_function(booster_x_exposures, booster_amount, power_H_x0, power_H_amount, target_x) + center_x;
	std::cout << "  First position of hammer x Global: " << hammer_x1 << std::endl;

	return 0;
}