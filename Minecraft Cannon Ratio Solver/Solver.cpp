#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>


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
        double f = (-1.0 / 8.0 + 1.0 / distance) * source.amount * exposure;
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
            f = -1.0 / 8.0 + 1.0 / distance;
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

    void printPos(std::string str = "") {
        std::cout << str << "x: " << x << "\ty: " << y << "\tz: " << z << std::endl;
    }
};
// Define static const members outside the class
const double Tnt::explosion_height = 0.061250001192092896;
const double Tnt::gravity = -0.04f;
const double Tnt::drag = 0.98f;

// Hash function for the vector
struct VectorHash {
    std::size_t operator()(const std::vector<int>& vec) const {
        std::size_t seed = 0;
        for (int num : vec) {
            seed ^= std::hash<int>()(num) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

bool generateNextPermutation(short int currentArray[], const short int startArray[], const short int finalArray[], size_t size) {
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

double fastff(double vel) {
    double start = 253.97999998182058;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    vel *= Tnt::drag;
    vel += Tnt::gravity;
    start += vel;
    return start;
}

std::mutex output_mutex;
double global_best_error = std::numeric_limits<double>::max(); // Start with the worst possible error

void process_range(int threadNum, const std::vector<int>& assigned_numbers) {
    std::cout << "Thread[" << threadNum << "] assigned numbers: ";
    unsigned long long int outputModNumber = 0;
    for (int num : assigned_numbers) {
        std::cout << num << " ";
    }
    std::cout << std::endl;

    for (int param : assigned_numbers) {
		// Reset/initialize has map for each thread
        std::cout << "Processing parameter: " << param << std::endl;

        Tnt swingTnt("-199974.5000002721 254.0199999809265 -333625.50999999046 0.9701997615808802 0.0 0.0", 75);
        Tnt swingTntRev("-199974.49000000954 248.0 -333625.49000000954", 1);
        swingTnt.explosion(swingTntRev, 1);

        Tnt powerTnt0("-199974.49755261914 254.0199999809265 -333620.5534764172 -0.00740155756734615 0.0 0.9177931383509221");
        Tnt powerTnt0Rev0("-199966.50999999046  254.0 -333620.11500000954", 1);
        powerTnt0.explosion(powerTnt0Rev0);

        std::vector<int> powerAdd;
        int add = -1;
        std::vector<double> xlocations, ylocations;
        int cap = 0;
        for (int i = param; i <= 1000; i++) {
            add++;
            Tnt powerTnt1 = powerTnt0; // reset
            swingTnt.amount = i;
            powerTnt1.swing(swingTnt);
            powerTnt1.freefall(1);
            powerTnt1.y = 253.98000010565525;

            double decimal = powerTnt1.x - std::floor(powerTnt1.x);
            if (decimal < 0.49f || decimal > 0.51f) {
                continue;
            }

            powerTnt1.z = -333592.49000000954; // guider coordinate
            powerAdd.push_back(add);
            powerTnt1.print("  swing + ff: [" + std::to_string(i) + "] (+" + std::to_string(add) + ") power: ");
            xlocations.push_back(powerTnt1.x);
            ylocations.push_back(powerTnt1.y);
            cap++;
            if (cap >= 15) { break; }
            add = 0;
        }

        std::vector<double> microExposures(15*13, 0.0);
        Tnt* microBoosters[15];
        for (int i = 0; i < 15; i++) {
            microBoosters[i] = new Tnt(xlocations[i], ylocations[i], -333592.49000000954, 0.0, 0.0, 0.0, 1);
        }

        Tnt projectile0("-199974.50999999046 253.97999998182058 -333592.49000000954 0.0 -0.03919999988675116 0.0");
        for (int i = 0; i < 15; i++) {
            Tnt projectile1i(projectile0.x, projectile0.y, projectile0.z, 0.0, 0.0, 0.0);
            microExposures[i * 13] = 0.0;
            for (int j = 1; j <= 12; j++) {
                projectile1i.explosion(*microBoosters[i]);
                microExposures[i * 13 + j] = projectile1i.v;
            }
        }
        Tnt yPower0("-199974.49000000954 0.0 -333592.49000000954 0.0 0.0 0.0");
        yPower0.y = 253.97998101488093;

        double errorYvel = std::numeric_limits<double>::max();
        double startVel255 = 0.0;

        for (int i = 0; i < 100000; i++) {
            yPower0.amount = i;
            Tnt projectile1 = projectile0;
            projectile1.explosion(yPower0);
            double errorCalculated = std::abs(projectile1.v - 0.31792389011258693);
            if (errorCalculated < errorYvel && projectile1.v > 0.31792389011258693) {
                errorYvel = errorCalculated;
                startVel255 = projectile1.v;
                std::cout << "  power: " << yPower0.amount << std::endl;
            }
        }
		std::this_thread::sleep_for(std::chrono::seconds(10)); // Add delay to stagger threads
        short int microAmountStart[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        short int microAmountCurrent[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        short int microAmountFinal[15] = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 };
        double finalVelocity255 = 0.0;
        double scoreY = 1.0;

        do {
            outputModNumber++;
			if (outputModNumber % 10000000000 == 0) {
				std::lock_guard<std::mutex> lock(output_mutex);
				std::cout << "Thread[" << threadNum << "] processed " << outputModNumber << " permutations." << std::endl;
				std::cout << "  Current vector searched: ";
				for (int i = 0; i < 15; i++) {
					std::cout << microAmountCurrent[i] << " ";
				}
                std::cout << std::endl;
			}
            finalVelocity255 = startVel255;
			for (int i = 0; i < 15; i++) {
				finalVelocity255 += microExposures[i * 13 + microAmountCurrent[i]];
			}
            //scoreY = fastff(finalVelocity255);
            if (finalVelocity255 < 0.3179238901125882) {
                // Check if all numbers to the right of the current-changing index are zero
                short int allZerosRight = 0;
                for (int j = 14; j >= 0; --j) {
                    if (microAmountCurrent[j] == 0) {
                        allZerosRight++;
                    }
                    else {
                        break;
                    }
                }
                microAmountCurrent[14 - allZerosRight] = microAmountFinal[14 - allZerosRight];
                continue; // Skip invalid scores but do not apply the finalization logic
            }
            if (finalVelocity255 == 0.3179238901125882) {
				std::cout << "Thread[" << threadNum << "] found a new best error: " << scoreY - 255.0  << std::endl;
                std::lock_guard<std::mutex> lock(output_mutex);
                global_best_error = scoreY - 255.0;
                std::cout << "  ";
                for (int i = 0; i < 15; i++) {
                    std::cout << microAmountCurrent[i] << " ";
                }
                std::cout << std::endl;
                std::cout << "  record255Y7: " << scoreY << std::endl;
                std::cout << "  Velocity255V0: " << finalVelocity255 << std::endl;
                std::cout << "  error255: " << global_best_error << std::endl;

                Tnt output(0.0, projectile0.y, 0.0, 0.0, finalVelocity255, 0.0);
                output.print("  tick 0: ");
                for (int i = 0; i < 7; i++) {
                    output.freefall(1);
                    output.print("  tick " + std::to_string(i) + ": ");
                }
                std::cout << "  final position: " << output.y << std::endl;
                std::cout << "  ratio input help: ";
                for (int i = 0; i < 15; i++) {
                    std::cout << powerAdd[i] << ":" << microAmountCurrent[i] << " ";
                }
                std::cout << std::endl;
            }
        } while (generateNextPermutation(microAmountCurrent, microAmountStart, microAmountFinal, 15));

        for (int i = 0; i < 15; i++) {
            delete microBoosters[i];
        }
    }

    std::cout << "Thread[" << threadNum << "] finished processing." << std::endl;
}

int main() {
    std::cout << std::setprecision(17);
    unsigned int thread_count = std::thread::hardware_concurrency();
    std::cout << "Available threads: " << thread_count << std::endl;

    const int start = 70, end = 210;
    std::vector<std::thread> threads;
    std::vector<std::vector<int>> thread_assignments(thread_count);

    int current_thread = 0;
    for (int i = start; i <= end; i += 10) {
        thread_assignments[current_thread].push_back(i);
        current_thread = (current_thread + 1) % thread_count;
    }

    for (unsigned int i = 0; i < thread_count; i++) {
        threads.emplace_back([i, &thread_assignments]() {
            std::this_thread::sleep_for(std::chrono::seconds(1) * i * 0.5); // Add delay to stagger threads
            process_range(i, thread_assignments[i]);
            });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    std::cout << "Final Best Error: " << global_best_error << std::endl;
    std::cout << "End of computation." << std::endl;
    while (true) {

    }
    return 0;
}