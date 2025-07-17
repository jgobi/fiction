#include "fiction_experiments.hpp"

#include <fiction/algorithms/simulation/sidb/operational_domain.hpp>      // operational domain computation algorithms
#include <fiction/algorithms/simulation/sidb/critical_temperature.hpp>
#include <fiction/algorithms/simulation/sidb/detect_bdl_pairs.hpp>
#include <fiction/algorithms/simulation/sidb/exhaustive_ground_state_simulation.hpp>
#include <fiction/algorithms/simulation/sidb/sidb_simulation_parameters.hpp>
#include <fiction/io/read_sqd_layout.hpp>
#include <fiction/io/write_operational_domain.hpp>                            // writer for operational domains
#include <fiction/io/print_layout.hpp>
#include <fiction/types.hpp>
#include <fiction/utils/truth_table_utils.hpp>

#include <fmt/format.h>

#include <array>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

using namespace fiction;

int main()  // NOLINT
{
    experiments::experiment<std::string, double, double, double, double, double, uint64_t, double, uint64_t, double, uint64_t, double, uint64_t, double> simulation_exp{
        "Benchmark",
        "Gate Name",
        "mu",
        "Critical Temperature [K] (old values)",
        "E_{g,err} [meV] (old values)",
        "Critical Temperature [K] (new values)",
        "E_{g,err} [meV] (new values)",
        "#Samples (CT)",
        "op. (CT)",
        "sim calls (CT)",
        "t in s (CT)",  // Contour Tracing
        "#Samples (FF)",
        "op. (FF)",
        "sim calls (FF)",
        "t in s (FF)",  // Flood Fill
    };

    static std::vector<std::pair<std::string, std::vector<tt>>> gates;

    static const std::string folder = fmt::format("{}qcastyle/gates/", EXPERIMENTS_PATH);
    for (const auto & entry : std::filesystem::directory_iterator(folder)) {
        std::string path = entry.path();
        if (path.find(".sqd") != std::string::npos) {
            std::string gate = path.substr(path.find_last_of("/\\") + 1);
            gate = gate.substr(0, gate.find_last_of('.'));
            std::cout << "Gate: " << gate << std::endl;
            const std::string gate_type = gate.substr(0, gate.find_first_of("_"));
            if (gate_type == "or") {
                gates.push_back(std::make_pair(gate, std::vector<tt>{create_or_tt()}));
            } else if (gate_type == "wire") {
                gates.push_back(std::make_pair(gate, std::vector<tt>{create_id_tt()}));
            } else if (gate_type == "not") {
                gates.push_back(std::make_pair(gate, std::vector<tt>{create_not_tt()}));
            } else if (gate_type == "wire2") {
                gates.push_back(std::make_pair(gate, create_double_wire_tt()));
            } else if (gate_type == "maj") {
                const uint32_t num_vars = 3;

                // Create a dynamic truth table with 3 variables
                kitty::dynamic_truth_table ttmaj(num_vars);
                // Binary string representing the truth table
                std::string binary_string = "00101011"; // MAJ for dealing with QCA style representation

                // Create the truth table from the binary string
                kitty::create_from_binary_string(ttmaj, binary_string);

                gates.push_back(std::make_pair(gate, std::vector<tt>{ttmaj}));
            } else {
                throw std::runtime_error(fmt::format("Unknown gate type: {}", gate_type));
            }
        }
    }

    const std::array<double, 1> mus = { -0.32 };

    for (const auto mu : mus ) {
        const sidb_simulation_parameters sim_params_new(3, mu, 4.1, 1.8);
        is_operational_params      isop_params_new{
                sim_params_new,
                sidb_simulation_engine::CLUSTERCOMPLETE
        };

        const sidb_simulation_parameters sim_params_old(3, mu, 5.6, 5);
        is_operational_params      isop_params_old{
                sim_params_old,
                sidb_simulation_engine::CLUSTERCOMPLETE
        };
        isop_params_old.input_bdl_iterator_params
            = isop_params_new.input_bdl_iterator_params
            // = bdl_input_iterator_params{ detect_bdl_wires_params {3, detect_bdl_pairs_params { 2.3, 2.9 }}}; // uncomment for Samuel thesis' wire
            // = bdl_input_iterator_params{ detect_bdl_wires_params {2, detect_bdl_pairs_params { 1.5, 1.9 }}}; // uncomment for 3x3 QCA-styled wires
            = bdl_input_iterator_params{};                                                                   // uncomment for 2z2 QCA-styled wires



        const critical_temperature_params ct_params_new{ isop_params_new };
        const critical_temperature_params ct_params_old{ isop_params_old };
        // operational domain parameters
        operational_domain_params op_domain_params{ isop_params_old };

        op_domain_params.operational_params.op_condition = is_operational_params::operational_condition::TOLERATE_KINKS;

        op_domain_params.sweep_dimensions         = {{sweep_parameter::EPSILON_R}, {sweep_parameter::LAMBDA_TF}};
        op_domain_params.sweep_dimensions[0].min  = 1.0;
        op_domain_params.sweep_dimensions[0].max  = 10.0;
        op_domain_params.sweep_dimensions[0].step = 0.05;
        op_domain_params.sweep_dimensions[1].min  = 1.0;
        op_domain_params.sweep_dimensions[1].max  = 10.0;
        op_domain_params.sweep_dimensions[1].step = 0.05;
        
        for (const auto& [gate, truth_table] : gates)
        {
            auto layout = read_sqd_layout<sidb_100_cell_clk_lyt_siqad>(fmt::format("{}/{}.sqd", folder, gate));

            critical_temperature_stats ct_stats_old{};
            const auto ct_old =
                critical_temperature_gate_based<sidb_100_cell_clk_lyt_siqad>(layout, truth_table, ct_params_old, &ct_stats_old);
            std::cout << "Critical temperature (old): " << ct_old << std::endl;
            
            critical_temperature_stats ct_stats_new{};
            const auto ct_new =
                critical_temperature_gate_based<sidb_100_cell_clk_lyt_siqad>(layout, truth_table, ct_params_new, &ct_stats_new);
            std::cout << "Critical temperature (new): " << ct_new << std::endl;

            // Write critical temperature results to CSV
            const std::string ct_csv_file = fmt::format("{}critical_temperature_{:.2f}_{}.csv", folder, mu, gate);
            std::ofstream ct_csv(ct_csv_file);
            if (ct_csv.is_open()) {
                ct_csv << "Critical Temperature (Old),Critical Temperature (New)\n";
                ct_csv << ct_old << "," << ct_new << "\n";
                ct_csv.close();
            } else {
                throw std::runtime_error(fmt::format("Failed to open file: {}", ct_csv_file));
            }
            
            // write operational domain parameters
            write_operational_domain_params write_op_domain_params{};
            write_op_domain_params.non_operational_tag = "0";
            write_op_domain_params.operational_tag     = "1";
            write_op_domain_params.writing_mode        = write_operational_domain_params::sample_writing_mode::ALL_SAMPLES;

            operational_domain_stats op_domain_stats_ct{};
            operational_domain_stats op_domain_stats_ff{};

            const auto op_domain_ct =
                operational_domain_contour_tracing(layout, truth_table, 100, op_domain_params, &op_domain_stats_ct);

            write_operational_domain(op_domain_ct,
                fmt::format("{}operational_domain_contour_tracing_siqad_{:.2f}_{}.csv", folder, mu, gate),
                write_op_domain_params);


            const auto op_domain_ff =
                operational_domain_flood_fill(layout, truth_table, 250, op_domain_params, &op_domain_stats_ff);

            write_operational_domain(op_domain_ff,
                fmt::format("{}operational_domain_flood_fill_siqad_{:.2f}_{}.csv", folder, mu, gate),
                write_op_domain_params);
            
            const auto operational_percentage_ct =
                static_cast<double>(op_domain_stats_ct.num_operational_parameter_combinations) /
                static_cast<double>(op_domain_stats_ct.num_evaluated_parameter_combinations);
            const auto operational_percentage_ff =
                static_cast<double>(op_domain_stats_ff.num_operational_parameter_combinations) /
                static_cast<double>(op_domain_stats_ff.num_evaluated_parameter_combinations);

            simulation_exp(gate, mu,
                ct_old, ct_stats_old.energy_between_ground_state_and_first_erroneous, 
                ct_new, ct_stats_new.energy_between_ground_state_and_first_erroneous,

                // Contour Tracing
                op_domain_stats_ct.num_evaluated_parameter_combinations, operational_percentage_ct,
                op_domain_stats_ct.num_simulator_invocations, mockturtle::to_seconds(op_domain_stats_ct.time_total),
                // Flood Fill
                op_domain_stats_ff.num_evaluated_parameter_combinations, operational_percentage_ff,
                op_domain_stats_ff.num_simulator_invocations, mockturtle::to_seconds(op_domain_stats_ff.time_total)
            );
            simulation_exp.save();
            std::cout<<std::endl;
        }
    }
    simulation_exp.table();
    return EXIT_SUCCESS;
}
