#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <vector>
#include <thread>
#include <mutex>

const int ratio_size = 3;
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
    static const std::vector<double> distances_dx_dy_dz(const Tnt& from, const Tnt& to) {
        
        double dx = to.x - from.x;
        double dy = to.y - from.y;
        double dz = to.z - from.z;
        std::vector<double> distances = { dx, dy, dz };
        return distances;
    }
    static const void print_distances(const Tnt& from, const Tnt& to, std::string name = "") {
        std::vector<double> distances = distances_dx_dy_dz(from, to);
        std::string space = "";
        if (name != "") {
            std::cout << name << std::endl;
            space = "\t";
        }
        std::cout << space + "dx: " << distances[0] << std::endl;
        std::cout << space + "dy: " << distances[1] << std::endl;
        std::cout << space + "dz: " << distances[2] << std::endl;
    }

    Tnt() = default;

    Tnt(const std::string& attributes, unsigned int explosivePower = 1.0) : amount(explosivePower) {
        std::istringstream stream(attributes);
        stream >> x >> y >> z >> u >> v >> w;
    };

    Tnt(double x, double y, double z, double u = 0.0, double v = 0.0, double w = 0.0, unsigned int explosivePower = 1.0) : x(x), y(y), z(z), u(u), v(v), w(w), amount(explosivePower) {};

    void explosion(Tnt source, double exposure = 1.0) {
        double dx = x - source.x;
        double dy = y - (source.y + explosion_height); // double dy = y - source.y; y-(source.y + explosion_height) 
        double dz = z - source.z;
        double distance = std::sqrt(dx * dx + dy * dy + dz * dz); // float

        if (distance > 8.0 || distance == 0.0) {
            return;
        }
        double f = (-1.0 / 8.0 + 1.0 / distance) * source.amount * exposure;

        dy = y - (source.y + explosion_height);
        u += f * dx;
        v += f * dy;
        w += f * dz;
    }

    void swing(Tnt source) {
        Tnt origin(source);
        std::vector<std::vector<double>> locations(source.amount);
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
            locations[i] = loc;
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
    void freefall_print(int ticks, std::string str = "") {
        print(str + " 0 ");
        for (int i = 0; i < ticks; i++) {
            v += Tnt::gravity;

            x += u;
            y += v;
            z += w;
            print(str + " " + std::to_string(i + 1) + " ");
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
const double Tnt::gravity = -0.04; // -0.04f
const double Tnt::drag = 0.98; //0.98f

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

bool generateNextPermutation(short int currentArray[], const short int startArray[], const short int finalArray[]) {
    size_t i = ratio_size;

    while (i > 0) {
        --i;
        if (currentArray[i] < finalArray[i]) {
            ++currentArray[i];
            for (size_t j = i + 1; j < ratio_size; ++j) {
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
double global_best_value = std::numeric_limits<double>::max(); // Start with the worst possible error

double min_height, max_height;
unsigned int final_hammer_amount;

void process_range(int threadNum, const std::vector<int>& assigned_powers) {
    // ---------------- config ----------------
    // Per-entity target heights (customize these)

    /*
    double target_y_sand = -62.036;
    double target_y_srh = -62.3311;
    double target_y_hr = -62.83;
    */

    // Fraction window for z
    const double frac_lo = 0.49, frac_hi = 0.51;
    const double distance_herv_threshold = 2.0;

    const double x_min_sand = -129092.5, x_max_sand = -129090.5;
    const double x_min_srh = -129102.5, x_max_srh = -129100.5;
    const double x_min_hr = -129112.5, x_max_hr = -129110.5;

    const double z0_sand = 67098.50999999046;
    const double z0_srh = 67098.50999999046;
    const double z0_hr = 67098.50999999046;

    const double Z_STEP = 1.0 / 16.0;
    const int Z_WINDOW_BLOCKS_SAND = 2;
    const int Z_WINDOW_BLOCKS_SRH = 2;
    const int Z_WINDOW_BLOCKS_HR = 2;

    const int    booster_amount = 40;
    const double x_tol = 1e-9;
    const int    max_bisect_iters = 100;

    const double power_y = -63.0; // "above power" means final y > -63.0
    const char* POWER_STR = "-129115.50999999046 -63.0 67096.5";

    const double SAND_X0 = -129092.49000000954, SAND_Y0 = -40.980000019073486;
    const double SRH_X0 = -129102.49000000954, SRH_Y0 = -40.980000019073486;
    const double HREV_X0 = -129112.49000000954, HREV_Y0 = -40.980000019073486;

    const double B_SAND_X0 = -129091.49000000954, B_Y0 = -40.980000019073486;
    const double B_SRH_X0 = -129100.64497978869;
    const double B_HR_X0 = -129110.67090365867;

    // Flight times
	const int FALL_TICKS = 11;
	const int SAND_FINAL_TICKS = 5;
	const int SRH_FINAL_TICKS = 6;
	const int HREV_FINAL_TICKS = 7;

    auto frac = [](double z) { return (z - std::floor(z)); };

    // -------- helpers that write to a generic ostream --------
    auto header = [](std::ostream& os, const std::string& title) {
        os << "\n========== " << title << " ==========\n";
        };
    auto kv = [](std::ostream& os, const std::string& k, double v) {
        os << "  " << std::left << std::setw(16) << k << ": " << std::setprecision(15) << v << "\n";
        };
    auto kvs = [](std::ostream& os, const std::string& k, const std::string& v) {
        os << "  " << std::left << std::setw(16) << k << ": " << v << "\n";
        };
    auto print_entity = [&](std::ostream& os, const char* name, const Tnt& e) {
        os << "  [" << name 
            << "] z=" << std::setprecision(15) << e.z
            << "  y=" << e.y
            << "  frac=" << frac(e.z)
            << "  block=" << (long)std::floor(e.z) << "\n";
        };
    auto print_paths_csv = [&](std::ostream& os, const std::string& label,
        const std::vector<std::pair<double, double>>& path) {
            // os << label << " (z,y)\n";
            for (size_t i = 0; i < path.size(); ++i) {
                os << path[i].first << "," << path[i].second << "\n";
            }
        };

    auto freefall_trace = [](Tnt& e, int ticks, std::vector<std::pair<double, double>>& out) {
        for (int i = 0; i < ticks; ++i) {
            e.v += Tnt::gravity;
            e.x += e.u;
            e.y += e.v;
            e.z += e.w;
            out.emplace_back(e.z, e.y);
            e.u *= Tnt::drag;
            e.v *= Tnt::drag;
            e.w *= Tnt::drag;
        }
        };

    // Simulate ONE entity end-to-end (booster -> FALL_TICKS ticks -> power -> final ticks)
    auto simulate_entity = [&](int entity_id, double booster_x, double z_snap, int n_pow, Tnt& out_entity) {
        Tnt power(POWER_STR);

        Tnt sand(SAND_X0, SAND_Y0, z_snap);
        Tnt srh(SRH_X0, SRH_Y0, z_snap);
        Tnt hrev(HREV_X0, HREV_Y0, z_snap);

        Tnt b_sand(B_SAND_X0, B_Y0, z_snap);
        Tnt b_srh(B_SRH_X0, B_Y0, z_snap);
        Tnt b_hr(B_HR_X0, B_Y0, z_snap);
        b_sand.amount = b_srh.amount = b_hr.amount = booster_amount;

        if (entity_id == 0) { b_sand.x = booster_x; sand.explosion(b_sand); }
        else if (entity_id == 1) { b_srh.x = booster_x; srh.explosion(b_srh); }
        else { b_hr.x = booster_x; hrev.explosion(b_hr); }

        // "Shot out" then coast
        sand.freefall(FALL_TICKS);
        srh.freefall(FALL_TICKS);
        hrev.freefall(FALL_TICKS);

        // Snap to power x and zero u before power explosion
        sand.x = -129115.50999999046; sand.u = 0.0;
        srh.x = -129115.50999999046; srh.u = 0.0;
        hrev.x = -129115.50999999046; hrev.u = 0.0;

        power.amount = n_pow;
        sand.explosion(power);
        srh.explosion(power);
        hrev.explosion(power);

        // Final ticks
        sand.freefall(SAND_FINAL_TICKS);
        srh.freefall(SRH_FINAL_TICKS);
        hrev.freefall(HREV_FINAL_TICKS);

        if (entity_id == 0) out_entity = sand;
        else if (entity_id == 1) out_entity = srh;
        else out_entity = hrev;
        };

    // Simulate ALL THREE with traces (booster -> FALL_TICKS -> power -> final ticks)
    auto simulate_full_chain_with_trace = [&](
        int n_pow,
        double x_sand, double z_sand,
        double x_srh, double z_srh,
        double x_hr, double z_hr,
        Tnt& sand_out, Tnt& srh_out, Tnt& hrev_out,
        Tnt& b_sand_out, Tnt& b_srh_out, Tnt& b_hr_out,
        std::vector<std::pair<double, double>>& sand_path,
        std::vector<std::pair<double, double>>& srh_path,
        std::vector<std::pair<double, double>>& hrev_path)
        {
            Tnt power(POWER_STR);

            Tnt sand(SAND_X0, SAND_Y0, z_sand);
            Tnt srh(SRH_X0, SRH_Y0, z_srh);
            Tnt hrev(HREV_X0, HREV_Y0, z_hr);

            Tnt b_sand(B_SAND_X0, B_Y0, z_sand);
            Tnt b_srh(B_SRH_X0, B_Y0, z_srh);
            Tnt b_hr(B_HR_X0, B_Y0, z_hr);
            b_sand.amount = b_srh.amount = b_hr.amount = booster_amount;

            b_sand.x = x_sand;
            b_srh.x = x_srh;
            b_hr.x = x_hr;

            // Booster explosions ("shot out")
            sand.explosion(b_sand);
            srh.explosion(b_srh);
            hrev.explosion(b_hr);

            // FALL_TICKS ticks coasting, capture traces
            freefall_trace(sand, FALL_TICKS, sand_path);
            freefall_trace(srh, FALL_TICKS, srh_path);
            freefall_trace(hrev, FALL_TICKS, hrev_path);

            // Snap to power line and zero horizontal u
            sand.x = -129115.50999999046; sand.u = 0.0;
            srh.x = -129115.50999999046; srh.u = 0.0;
            hrev.x = -129115.50999999046; hrev.u = 0.0;

            // Apply power
            power.amount = n_pow;
            sand.explosion(power);
            srh.explosion(power);
            hrev.explosion(power);

            // Final ticks and traces
            freefall_trace(sand, SAND_FINAL_TICKS, sand_path);
            freefall_trace(srh, SRH_FINAL_TICKS, srh_path);
            freefall_trace(hrev, HREV_FINAL_TICKS, hrev_path);

            sand_out = sand; srh_out = srh; hrev_out = hrev;
            b_sand_out = b_sand; b_srh_out = b_srh; b_hr_out = b_hr;
        };

    // Final y helper for root finding against a target
    auto final_y_entity = [&](int entity_id, double booster_x, double z_snap, int n_pow) -> double {
        Tnt out_e;
        simulate_entity(entity_id, booster_x, z_snap, n_pow, out_e);
        return out_e.y;
        };

    // Objective for bisection: final_y - target_y
    auto f_eval = [&](int entity_id, double x, double z_snap, int n_pow, double target_y) -> double {
        return final_y_entity(entity_id, x, z_snap, n_pow) - target_y;
        };

    // Find an x that yields final y below target_y; if entire interval is below, pick the endpoint
    // whose final y is closest to target from below (near-but-under behavior).
    auto solve_booster_x_below = [&](
        int entity_id,
        double x_min, double x_max,
        double z_snap,
        int n_pow,
        double target_y,
        double& x_star, double& y_star,
        std::string& reason) -> bool
        {
            double f_lo = f_eval(entity_id, x_min, z_snap, n_pow, target_y);
            double f_hi = f_eval(entity_id, x_max, z_snap, n_pow, target_y);

            if (f_lo >= 0.0 && f_hi >= 0.0) {
                reason = "no x gives final y below target";
                return false;
            }
            if (f_lo < 0.0 && f_hi < 0.0) {
                double y_lo = final_y_entity(entity_id, x_min, z_snap, n_pow);
                double y_hi = final_y_entity(entity_id, x_max, z_snap, n_pow);
                if (y_lo > y_hi) { x_star = x_min; y_star = y_lo; }
                else { x_star = x_max; y_star = y_hi; }
                reason = "entire x-range below target; chose endpoint closest to target (near but under)";
                return true;
            }

            double lo = x_min, hi = x_max;
            double flo = f_lo, fhi = f_hi;
            if (!(flo < 0.0 && fhi >= 0.0)) {
                std::swap(lo, hi);
                std::swap(flo, fhi);
            }

            int it = 0;
            while ((hi - lo) > x_tol && it < max_bisect_iters) {
                double mid = 0.5 * (lo + hi);
                double fmid = f_eval(entity_id, mid, z_snap, n_pow, target_y);
                if (fmid < 0.0) { lo = mid; flo = fmid; }
                else { hi = mid; fhi = fmid; }
                ++it;
            }

            x_star = lo;
            y_star = final_y_entity(entity_id, x_star, z_snap, n_pow);
            reason.clear();
            return true;
        };

    auto meets_constraints = [&](const Tnt& e, double target_y) -> bool {
        double fz = frac(e.z);
        return (e.y < target_y) && (e.y > power_y) && (fz >= frac_lo) && (fz <= frac_hi);
        };

    {
        std::ostringstream oss;
        oss << "Thread [" << threadNum << "] assigned powers: ";
        for (int np : assigned_powers) {
            oss << np << " ";
        }
        oss << "\n";

        std::lock_guard<std::mutex> lk(output_mutex);
        std::cout << oss.str();
    }
	// Pause a moment to stagger thread starts (helps with console output readability)
    std::this_thread::sleep_for(std::chrono::seconds(1));

    // =================== MAIN SEARCH ===================
    for (int n_pow : assigned_powers) {
        {
            std::ostringstream oss;
            oss << "Thread [" << threadNum << "] processing n_pow = " << n_pow << "\n";
            std::lock_guard<std::mutex> lk(output_mutex);
            std::cout << oss.str();
        }

        short int start_ratio[ratio_size] = { -4, -4, -4 };
        short int current_ratio[ratio_size] = { -4, -4, -4 };
        short int final_ratio[ratio_size] = { 4, 4, 4 };
        do {
            double target_y_sand = -62.1 + current_ratio[0] * 0.002;
            double target_y_srh = -62.3 + current_ratio[1] * 0.002;
            double target_y_hr = -62.8 + current_ratio[2] * 0.002;


            // Stage 1: SAND
            for (int iz_s = -16 * Z_WINDOW_BLOCKS_SAND; iz_s <= 16 * Z_WINDOW_BLOCKS_SAND; ++iz_s) {
                double z_sand = std::round((z0_sand + iz_s * Z_STEP) * 16.0) / 16.0;
                double x_sand, y_sand; std::string why_sand;
                if (!solve_booster_x_below(0, x_min_sand, x_max_sand, z_sand, n_pow,
                    target_y_sand, x_sand, y_sand, why_sand))
                    continue;

                Tnt sand_fin;
                simulate_entity(0, x_sand, z_sand, n_pow, sand_fin);
                if (!meets_constraints(sand_fin, target_y_sand)) continue;

                // Stage 2: SR_H (must share z block with SAND)
                for (int iz_r = -16 * Z_WINDOW_BLOCKS_SRH; iz_r <= 16 * Z_WINDOW_BLOCKS_SRH; ++iz_r) {
                    double z_srh = std::round((z0_srh + iz_r * Z_STEP) * 16.0) / 16.0;
                    double x_srh, y_srh; std::string why_srh;
                    if (!solve_booster_x_below(1, x_min_srh, x_max_srh, z_srh, n_pow,
                        target_y_srh, x_srh, y_srh, why_srh))
                        continue;

                    Tnt srh_fin;
                    simulate_entity(1, x_srh, z_srh, n_pow, srh_fin);
                    if (!meets_constraints(srh_fin, target_y_srh)) continue;

                    long bs = (long)std::floor(sand_fin.z);
                    long br = (long)std::floor(srh_fin.z);
                    if (bs != br) continue;

                    // Stage 3: HR (final)
                    for (int iz_h = -16 * Z_WINDOW_BLOCKS_HR; iz_h <= 16 * Z_WINDOW_BLOCKS_HR; ++iz_h) {
                        double z_hr = std::round((z0_hr + iz_h * Z_STEP) * 16.0) / 16.0;
                        double x_hr, y_hr; std::string why_hr;
                        if (!solve_booster_x_below(2, x_min_hr, x_max_hr, z_hr, n_pow,
                            target_y_hr, x_hr, y_hr, why_hr))
                            continue;

                        Tnt hr_fin;
                        simulate_entity(2, x_hr, z_hr, n_pow, hr_fin);
                        if (Tnt::distance(hr_fin, srh_fin) > distance_herv_threshold && (hr_fin.y + Tnt::explosion_height) < srh_fin.y) continue;

                        // Full chain found — gather traces and print ONCE atomically
                        Tnt sand_f, srh_f, hrev_f, b_sand_f, b_srh_f, b_hr_f;
                        std::vector<std::pair<double, double>> sand_path, srh_path, hrev_path;
                        sand_path.reserve(FALL_TICKS + SAND_FINAL_TICKS);
                        srh_path.reserve(FALL_TICKS + SRH_FINAL_TICKS);
                        hrev_path.reserve(FALL_TICKS + HREV_FINAL_TICKS);

                        simulate_full_chain_with_trace(
                            n_pow,
                            x_sand, z_sand,
                            x_srh, z_srh,
                            x_hr, z_hr,
                            sand_f, srh_f, hrev_f,
                            b_sand_f, b_srh_f, b_hr_f,
                            sand_path, srh_path, hrev_path
                        );

                        std::ostringstream oss;
                        header(oss, "H-HIT (FULL CHAIN + PATHS)");
                        kv(oss, "n_pow", n_pow);
                        kv(oss, "S.z_snap", z_sand); kv(oss, "S.booster_x", x_sand);
                        kv(oss, "R.z_snap", z_srh);  kv(oss, "R.booster_x", x_srh);
                        kv(oss, "H.z_snap", z_hr);   kv(oss, "H.booster_x", x_hr);

                        oss << "\n-- Final Projectile States --\n";
                        print_entity(oss, "S", sand_f);
                        print_entity(oss, "SR_H", srh_f);
                        print_entity(oss, "HR", hrev_f);

                        oss << "\n-- Booster Projectile --\n";
                        print_entity(oss, "B_S", b_sand_f);
                        print_entity(oss, "B_SRH", b_srh_f);
                        print_entity(oss, "B_HR", b_hr_f);

                        oss << "\n-- Trajectory Data --\n";
                        oss << "  Distance Between Projectiles:"
                            << "\n    (S-R): " << Tnt::distance(sand_f, srh_f)
                            << "\n    (S-H): " << Tnt::distance(sand_f, hrev_f)
                            << "\n    (R-H): " << Tnt::distance(srh_f, hrev_f);
                        Tnt tempPower(POWER_STR);
                        oss << "\n  Range: " << sand_f.z - tempPower.z << " meters = " << (sand_f.z - tempPower.z) / 16.0 << " m16\n";

                        oss << "\n-- Flight Times --\n";
                        oss << "  Freefall: " << FALL_TICKS << "\n";
                        oss << "  S final: " << SAND_FINAL_TICKS << "\n";
                        oss << "  SR_H final: " << SRH_FINAL_TICKS << "\n";
                        oss << "  HR final: " << HREV_FINAL_TICKS << "\n";

                        oss << "\n-- Flight Paths (z,y) --\n";
                        print_paths_csv(oss, "S_path", sand_path);
                        print_paths_csv(oss, "SR_H_path", srh_path);
                        print_paths_csv(oss, "HR_path", hrev_path);
                        oss << "-----------------------------------------\n";


                        // want to find how much rev tnt is needed for sand rev and hammer rev to push sand and hammer vertically up to same block (within y.0 and y.02)
                        Tnt sandReved = sand_f;
                        Tnt sandRev = srh_f;
                        Tnt hammerReved = srh_f;
                        Tnt hammerRev = hrev_f;

						bool found = false;

                        for (int srev_tnt = 0; srev_tnt <= 500; srev_tnt++) {
                            for (int hrev_tnt = 0; hrev_tnt <= 500; hrev_tnt++) {
                                sandReved = sand_f;
                                sandRev = srh_f;
                                hammerReved = srh_f;
                                hammerRev = hrev_f;

                                sandRev.amount = srev_tnt;
                                hammerRev.amount = hrev_tnt;

                                sandReved.explosion(sandRev);
                                hammerReved.explosion(hammerRev);

                                sandReved.freefall(1);
                                hammerReved.freefall(1);

                                double decimalSand = hammerReved.y - std::floor(hammerReved.y);
                                int blockSand = (int)std::floor(hammerReved.y);
                                double decimalHammer = sandReved.y - std::floor(sandReved.y);
                                int blockHammer = (int)std::floor(sandReved.y);
								double sand_dz = sandReved.z - sand_f.z;
								double hammer_dz = hammerReved.z - srh_f.z;

                                if (blockSand == blockHammer && decimalSand > 0.0 && decimalSand < 0.02 && decimalHammer > 0.0 && decimalHammer < 0.02 && sand_dz > 0 && hammer_dz > 0) {// && decimalSand > 0.0 && decimalSand < 0.02 && decimalHammer > 0.0 && decimalHammer < 0.02
                                    found = true;
                                    oss << "  S REV EXP: " << srev_tnt << " (final z,y: " << std::setprecision(15) << sandReved.z << ", " << sandReved.y << ") dz: " << sand_dz << "\n";
                                    oss << "  HR REV EXP: " << hrev_tnt << " (final z,y: " << std::setprecision(15) << hammerReved.z << ", " << hammerReved.y << ") dz: " << hammer_dz << "\n";
                                    oss << "\n";
                                }
                            }
                        }

                        if (found)
                        {
                            std::lock_guard<std::mutex> lk(output_mutex);
                            std::cout << oss.str();
                        }
                        else {
                            oss.clear();
                        }
                    }
                }
            }
		} while (generateNextPermutation(current_ratio, start_ratio, final_ratio));
    }
}

int main() {
    std::cout << std::setprecision(17);

    // Configure the global power range here
    const int n_pow_min = 40;
    const int n_pow_max = 200;

    unsigned int thread_count = std::max(1u, std::thread::hardware_concurrency() - 1u);
    std::cout << "Available threads: " << thread_count << std::endl;

    const int start = n_pow_min, end = n_pow_max;
    std::vector<std::thread> threads;
    std::vector<std::vector<int>> thread_assignments(thread_count);

    int current_thread = 0;
    for (int i = start; i <= end; i++) {
        thread_assignments[current_thread].push_back(i);
        current_thread = (current_thread + 1) % thread_count;
    }

    for (unsigned int i = 0; i < thread_count; i++) {
        threads.emplace_back([i, &thread_assignments]() {
            std::this_thread::sleep_for(std::chrono::seconds(1) * i * 0.05); // Add delay to stagger threads
            process_range(i, thread_assignments[i]);
            });
    }

    for (auto& thread : threads) {
        thread.join();
    }

    // Print results
    std::cout << "Final Best Value: " << global_best_value << std::endl;
    std::cout << "End of computation." << std::endl;
    while (true) {

    }
    return 0;
}