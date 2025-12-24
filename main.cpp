#include <iostream>
#include <vector>
#include <string>
#include <matplotlibcpp.h>
#include <Python.h>
#include "solver.hpp"
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <cmath>
#include <chrono>


using namespace std;
namespace plt = matplotlibcpp;


int main(int argc, char* argv[]){
    if (argc < 2){
        cout << "Usage: " << argv[0] << " <netlist_file>\n";
        return 1;
    }
    string netlist_file = argv[1];
    circuit ckt;
    vector<analysis> analyses;
    parseNetlistFile(netlist_file, ckt, analyses);
    if (analyses.empty()) {
        cerr << "No analysis specified in the netlist.\n";
        return 1;
    }
    for (auto& analysis : analyses){
        solver sol(ckt, analysis);
        if (analysis.type == "DC"){
            //询问DC分析方法
            cout << "Select DC analysis method:\n";
            cout << "1. Gauss Elimination\n";
            cout << "2. LU Decomposition\n";
            cout << "3. Manual LU\n";
            cout << "4. Gauss-Jacobi Iteration\n";
            cout << "5. Gauss-Seidel Iteration\n";
            int method_choice;
            cin >> method_choice;
            switch (method_choice){
                case 1:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_ELIMINATION);
                    break;
                case 2:
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
                    break;
                case 3:
                    sol.setLinearSolverMethod(LinearSolverMethod::MANUAL_LU);
                    break;
                case 4:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_JACOBI);
                    break;
                case 5:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_SEIDEL);
                    break;
                default:
                    cout << "Invalid choice, using LU Decomposition by default.\n";
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
            }
            //解DC分析
            auto dc_start = std::chrono::steady_clock::now();
            sol.DC_solve_ramp();
            auto dc_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> dc_elapsed = dc_end - dc_start;
            cout << "DC solve time (s): " << dc_elapsed.count() << "\n";
            //打印要观察的节点电压
            for (int node_id : sol.get_plot_node_ids()){
                sol.print_node_voltage(node_id);
            }
            //打印要观察的支路电流
            for (int node_id : sol.get_plot_current_ids()){
                sol.print_branch_current(node_id);
            }
        }
        else if (analysis.type == "TRAN"){
            double tstep = analysis.parameters.count("tstep") ? analysis.parameters["tstep"] : 1e-9;
            double tstop = analysis.parameters.count("tstop") ? analysis.parameters["tstop"] : 1e-6;
            //选择瞬态分析方法
            cout << "Select Transient analysis method:\n";
            cout << "1. Forward Euler\n";
            cout << "2. Backward Euler\n";
            cout << "3. Trapezoidal\n";
            int tran_method_choice;
            cin >> tran_method_choice;
            switch (tran_method_choice){
                case 1:
                    sol.setTransientMethod(TransientMethod::FORWARD_EULER);
                    break;
                case 2:
                    sol.setTransientMethod(TransientMethod::BACKWARD_EULER);
                    break;
                case 3:
                    sol.setTransientMethod(TransientMethod::TRAPEZOIDAL);
                    break;
                default:
                    cout << "Invalid choice, using Trapezoidal method by default.\n";
                    sol.setTransientMethod(TransientMethod::TRAPEZOIDAL);
            }

            //选择零初值还是直流点作为初值
            cout << "Select initial condition for Transient analysis:\n";
            cout << "1. Zero Initial Condition\n";
            cout << "2. DC Operating Point\n";
            int ic_choice;
            cin >> ic_choice;
            //询问DC分析方法
            cout << "Select DC analysis method:\n";
            cout << "1. Gauss Elimination\n";
            cout << "2. LU Decomposition\n";
            cout << "3. Manual LU\n";
            cout << "4. Gauss-Jacobi Iteration\n";
            cout << "5. Gauss-Seidel Iteration\n";
            int method_choice;
            cin >> method_choice;
            switch (method_choice){
                case 1:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_ELIMINATION);
                    break;
                case 2:
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
                    break;
                case 3:
                    sol.setLinearSolverMethod(LinearSolverMethod::MANUAL_LU);
                    break;
                case 4:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_JACOBI);
                    break;
                case 5:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_SEIDEL);
                default:
                    cout << "Invalid choice, using LU Decomposition by default.\n";
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
            }

            //sol.TRAN_solve(tstop, tstep,ic_choice);
            //sol.TRAN_solve();
            //开始计时
            auto start_time = std::chrono::high_resolution_clock::now();
            sol.TRAN_solve(tstop, tstep, ic_choice);
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            cout << "Transient analysis completed in " << elapsed.count() << " seconds.\n";
            // //打印要观察的节点电压时间序列
            // for (const auto& tran_plot_pair : sol.get_tran_plot_data()){
            //     int node_id = tran_plot_pair.first;
            //     const auto& time_voltage_series = tran_plot_pair.second;
            //     cout << "Transient data for Node ID " << node_id << ":\n";
            //     for (const auto& tv : time_voltage_series){
            //         cout << "Time: " << tv.first << " s, Voltage: " << tv.second << " V\n";
            //     }
            // }

            //打印瞬态结果到文件
            sol.print_tran_results();

            // 12.14 绘制要观察的节点电压波形和支路电流波形，要求画在一个窗口里，电流电压坐标轴独立
            plt::figure(1);
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
                int node_id = tran_plot_pair.first; 
                const auto& time_voltage_series = tran_plot_pair.second; 
                vector<double> times, values; 
                for (const auto& tv : time_voltage_series){ 
                    times.push_back(tv.first); 
                    values.push_back(tv.second); 
                } 
                if (node_id >= 0 && node_id < (int)ckt.node_list.size()){
                    plt::plot(times, values, {{"label", "Node " + to_string(node_id) + "(" + ckt.node_list[node_id] + ") Voltage"}}); 
                }
                else {
                    continue;
                }
            } 
            plt::legend(); 
            plt::xlabel("Time (s)"); 
            plt::ylabel("Value (V)"); 
            plt::title("Transient Analysis Node Voltages"); 
            plt::grid(true); 
            plt::show();

            plt::figure(2);
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
                int node_id = tran_plot_pair.first; 
                const auto& time_voltage_series = tran_plot_pair.second; 
                vector<double> times, values; 
                for (const auto& tv : time_voltage_series){ 
                    times.push_back(tv.first); 
                    values.push_back(tv.second); 
                } 
                if (node_id >= 0 && node_id < (int)ckt.node_list.size()){
                    continue;
                }
                else {
                    std::cout << "Plotting current for branch current index: " << -node_id << ckt.sources[-node_id-1].name << "\n";
                    plt::plot(times, values, {{"label", "I(" + ckt.sources[-node_id-1].name + ") "}}); 
                }
            } 
            plt::legend(); 
            plt::xlabel("Time (s)"); 
            plt::ylabel("Current (A)"); 
            plt::title("Transient Analysis Branch Currents"); 
            plt::grid(true); 
            plt::show();



            // //绘制要观察的节点电压波形，要求画在一个窗口里，并且坐标轴独立
            // plt::figure();
            // for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
            //     int node_id = tran_plot_pair.first; 
            //     const auto& time_voltage_series = tran_plot_pair.second; 
            //     vector<double> times, voltages; 
            //     for (const auto& tv : time_voltage_series){ 
            //         times.push_back(tv.first); 
            //         voltages.push_back(tv.second); 
            //     } 
            //     plt::plot(times, voltages, {{"label", "Node " + to_string(node_id) + "(" + ckt.node_list[node_id] + ")"}}); 
            // } 
            // plt::legend(); 
            // plt::xlabel("Time (s)"); 
            // plt::ylabel("Voltage (V)"); 
            // plt::title("Transient Analysis Node Voltages"); 
            // plt::grid(true); 
            // plt::show();
        }
        else if (analysis.type == "HB"){
            //询问使用0初值还是瞬态结果作为初值
            cout << "Select initial condition for Harmonic Balance analysis:\n";
            cout << "1. Zero Initial Condition\n";
            cout << "2. Use Transient Analysis Result\n";
            int ic_choice_hb;
            cin >> ic_choice_hb;
            //选择HB线性方程求解器
            cout << "Select HB linear solver:\n";
            cout << "1. Eigen LU (library)\n";
            cout << "2. Manual LU\n";
            int hb_solver_choice;
            cin >> hb_solver_choice;
            switch (hb_solver_choice){
                case 1:
                    sol.setHBLinearSolverMethod(HBLinearSolverMethod::EIGEN_LU);
                    break;
                case 2:
                    sol.setHBLinearSolverMethod(HBLinearSolverMethod::MANUAL_LU);
                    break;
                default:
                    cout << "Invalid choice, using Eigen LU by default.\n";
                    sol.setHBLinearSolverMethod(HBLinearSolverMethod::EIGEN_LU);
            }
            double freq = analysis.parameters.count("freq") ? analysis.parameters["freq"] : 1e3;
            double harm = analysis.parameters.count("harm") ? analysis.parameters["harm"] : 5;
            sol.PSS_solve_harmonic_balance(analysis,ic_choice_hb);
            //输出打印时域结果
            sol.print_hb_time_domain_results();
            cout << "\n";
            //输出打印频域结果
            sol.print_hb_frequency_domain_results();
            //输出打印要观察的节点电压
            for (int node_id : sol.get_plot_node_ids()){
                sol.print_node_voltage(node_id);
            }
            //绘制节点时域波形
            //sol.plot_hb_time_domain_results();

            // 我超级想知道为什么这段代码在solver.cpp里画不出来图。我知道了
            // 时间采样点数
            int N = 2 * sol.hb_params.num_harmonics + 1;
            double T = 1.0 / (sol.hb_params.fundamental_omega / (2.0 * M_PI)); //周期

            // 构造时间轴
            std::vector<double> time_points;
            for (int n = 0; n < N; ++n) {
                time_points.push_back(n * T / N);
            }

            // 创建一个图
            plt::figure();

            for (int node_id : sol.get_plot_node_ids()) {
                std::string name = "NODE";
                if (node_id >= 0 && node_id < (int)ckt.node_list.size())
                    name = ckt.node_list[node_id];
                // 提取节点电压
                std::vector<double> voltages;
                for (int n = 0; n < N; ++n) {
                    voltages.push_back(sol.hb_xt((node_id - 1) + n * sol.base_size).real());
                    //std::cout << time_points[n] << "\t" << voltages[n] << "\n";
                }
                // 绘制节点电压波形
                plt::plot(time_points, voltages, {{"label", "Node " + name}});
                //std::cout << "Plotted Node " << name << " voltage.\n";
            }
            plt::legend();
            plt::title("HB Time-Domain Node Voltages");
            plt::xlabel("Time (s)");
            plt::ylabel("Voltage (V)");
            plt::grid(true);
            plt::show();

            // 创建电流图
            plt::figure();

            for (int node_id : sol.get_plot_current_ids()) {
                std::string name = "NODE";
                if (node_id >= 0 && node_id < (int)ckt.sources.size())
                    name = ckt.sources[node_id].name;
                // 提取节点电流
                std::vector<double> currents;
                for (int n = 0; n < N; ++n) {
                    int branch_index = node_id + (ckt.node_list.size()-1);
                    currents.push_back(sol.hb_xt(branch_index + n * sol.base_size).real());
                    //std::cout << time_points[n] << "\t" << currents[n] << "\n";
                }
                // 绘制节点电流波形
                plt::plot(time_points, currents, {{"label", "I(" + name + ")"}});
                //std::cout << "Plotted Node " << name << " voltage.\n";
            }
            plt::legend();
            plt::title("HB Time-Domain Node Currents");
            plt::xlabel("Time (s)");
            plt::ylabel("Current (A)");
            plt::grid(true);
            plt::show();

            std::cout << "Plotted HB time-domain voltages for all nodes.\n";

        }
        else if (analysis.type == "SHOOTING"){
            double freq;
            if (analysis.parameters.count("freq")) {
                // 用户显式指定：保持原行为（你也可以改成“把它也纳入 gcd”，看你需求）
                freq = analysis.parameters["freq"];
            } else {
                std::vector<double> sin_freqs;
                sin_freqs.reserve(ckt.sources.size());

                for (const auto& source : ckt.sources) {
                    if (source.type == "V_SIN" || source.type == "I_SIN") {
                        auto it = source.parameters.find("FREQ");
                        if (it != source.parameters.end()) {
                            sin_freqs.push_back(it->second);
                        }
                    }
                }

                double f0 = 0.0;
                if (!gcd_fundamental_freq(sin_freqs, f0)) {
                    cerr << "No frequency specified for SHOOTING analysis and no valid sinusoidal source freq found.\n";
                    continue;
                }

                freq = f0;
                cout << "Can't find freq in input, initialize FREQ:" << freq << "\n";
            }

            //输入采样点数
            cout << "Enter number of samples for SHOOTING analysis (e.g., 1000): ";
            double N_sample;
            cin >> N_sample;

            double period_T = 1.0 / freq;
            double tstep = period_T / N_sample; //默认时间步长为周期的1/1000
            cout << "Select DC analysis method:\n";
            cout << "1. Gauss Elimination\n";
            cout << "2. LU Decomposition\n";
            cout << "3. Manual LU\n";
            cout << "4. Gauss-Jacobi Iteration\n";
            cout << "5. Gauss-Seidel Iteration\n";
            int dc_method_choice;
            cin >> dc_method_choice;
            switch (dc_method_choice){
                case 1:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_ELIMINATION);
                    break;
                case 2:
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
                    break;
                case 3:
                    sol.setLinearSolverMethod(LinearSolverMethod::MANUAL_LU);
                    break;
                case 4:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_JACOBI);
                    break;
                case 5:
                    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_SEIDEL);
                    break;
                default:
                    cout << "Invalid choice, using LU Decomposition by default.\n";
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
            }

            cout << "Enter pre-run cycles for SHOOTING (0 to skip): ";
            int pre_run_cycles;
            cin >> pre_run_cycles;
            cout << "Select the way of shooting method:\n";
            cout << "1. TR\n";
            cout << "2. BE\n";
            cout << "3. BE_sensitivity\n";
            cout << "4. TR_sensitivity\n";
            int method_choice;
            cin >> method_choice;
            auto shooting_start = std::chrono::steady_clock::now();
            switch (method_choice){
                case 1:
                    sol.PSS_solve_shooting_trapezoidal(period_T, tstep, 100, 1e-9, pre_run_cycles);
                    break;
                case 2:
                    sol.PSS_solve_shooting_backward_euler(period_T, tstep, 100, 1e-9, pre_run_cycles);
                    break;
                case 3:
                    sol.PSS_solve_shooting_backward_euler_sensitivity(period_T, tstep, 100, 1e-9, pre_run_cycles);
                    break;
                case 4:
                    sol.PSS_solve_shooting_trapezoidal_sensitivity(period_T, tstep, 100, 1e-9, pre_run_cycles);
                    break;
                default:
                    cout << "Invalid choice, using trapezoidal by default.\n";
                    sol.PSS_solve_shooting_trapezoidal(period_T, tstep, 100, 1e-9, pre_run_cycles);
            }
            auto shooting_end = std::chrono::steady_clock::now();
            std::chrono::duration<double> shooting_elapsed = shooting_end - shooting_start;
            cout << "Shooting solve time (s): " << shooting_elapsed.count() << "\n";

            // // Debug: 输出要plot的节点ID
            // cout << "Nodes to plot:\n";
            // for (int node_id : sol.get_plot_node_ids()){
            //     cout << "Node ID: " << node_id << "\n";
            // }
            // system("pause");

            // //打印要观察的节点电压
            // for (int node_id : sol.get_plot_node_ids()){
            //     sol.print_node_voltage(node_id);
            // }

            plt::figure(1);
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
                int node_id = tran_plot_pair.first; 
                const auto& time_voltage_series = tran_plot_pair.second; 
                vector<double> times, values; 
                for (const auto& tv : time_voltage_series){ 
                    times.push_back(tv.first); 
                    values.push_back(tv.second); 
                } 
                if (node_id >= 0 && node_id < (int)ckt.node_list.size()){
                    plt::plot(times, values, {{"label", "Node " + to_string(node_id) + "(" + ckt.node_list[node_id] + ") Voltage"}}); 
                }
                else {
                    continue;
                }
            } 
            plt::legend(); 
            plt::xlabel("Time (s)"); 
            plt::ylabel("Value (V)"); 
            plt::title("Shooting Analysis Node Voltages"); 
            plt::grid(true); 
            plt::show();

            plt::figure(2);
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
                int node_id = tran_plot_pair.first; 
                const auto& time_voltage_series = tran_plot_pair.second; 
                vector<double> times, values; 
                for (const auto& tv : time_voltage_series){ 
                    times.push_back(tv.first); 
                    values.push_back(tv.second); 
                } 
                if (node_id >= 0 && node_id < (int)ckt.node_list.size()){
                    continue;
                }
                else {
                    std::cout << "Plotting current for branch current index: " << -node_id << ckt.sources[-node_id-1].name << "\n";
                    plt::plot(times, values, {{"label", "I(" + ckt.sources[-node_id-1].name + ") "}}); 
                }
            } 
            plt::legend(); 
            plt::xlabel("Time (s)"); 
            plt::ylabel("Current (A)"); 
            plt::title("Shooting Analysis Branch Currents"); 
            plt::grid(true); 
            plt::show();

            // plt::figure();
            // for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
            //     int node_id = tran_plot_pair.first; 
            //     const auto& time_voltage_series = tran_plot_pair.second; 
            //     vector<double> times, voltages; 
            //     for (const auto& tv : time_voltage_series){ 
            //         times.push_back(tv.first); 
            //         voltages.push_back(tv.second); 
            //     } 
            //     plt::plot(times, voltages, {{"label", "Node " + to_string(node_id) + "(" + ckt.node_list[node_id] + ")"}}); 
            // } 
            // plt::legend(); 
            // plt::xlabel("Time (s)"); 
            // plt::ylabel("Voltage (V)"); 
            // plt::title("Transient Analysis Node Voltages"); 
            // plt::grid(true); 
            // plt::show();
        }
        else {
            cerr << "Unknown analysis type: " << analysis.type << "\n";
        }

    }
    return 0;
}
