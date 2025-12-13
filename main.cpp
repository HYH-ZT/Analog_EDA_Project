#include <iostream>
#include <vector>
#include <string>
#include <matplotlibcpp.h>
#include <Python.h>
#include "solver.hpp"
#include "circuit.hpp"
#include "parse_netlist.hpp"

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
                default:
                    cout << "Invalid choice, using LU Decomposition by default.\n";
                    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
            }
            //解DC分析
            sol.DC_solve();
            //打印要观察的节点电压
            for (int node_id : sol.get_plot_node_ids()){
                sol.print_node_voltage(node_id);
            }
        }
        else if (analysis.type == "TRAN"){
            double tstep = analysis.parameters.count("step") ? analysis.parameters["step"] : 1e-9;
            double tstop = analysis.parameters.count("stop") ? analysis.parameters["stop"] : 1e-6;
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
            sol.TRAN_solve(tstop, tstep);
            //sol.TRAN_solve();
            cout << "Transient analysis completed.\n";
            //打印要观察的节点电压时间序列
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){
                int node_id = tran_plot_pair.first;
                const auto& time_voltage_series = tran_plot_pair.second;
                cout << "Transient data for Node ID " << node_id << ":\n";
                for (const auto& tv : time_voltage_series){
                    cout << "Time: " << tv.first << " s, Voltage: " << tv.second << " V\n";
                }
            }
            //绘制要观察的节点电压波形，要求画在一个窗口里，并且坐标轴独立
            plt::figure();
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
                int node_id = tran_plot_pair.first; 
                const auto& time_voltage_series = tran_plot_pair.second; 
                vector<double> times, voltages; 
                for (const auto& tv : time_voltage_series){ 
                    times.push_back(tv.first); 
                    voltages.push_back(tv.second); 
                } 
                plt::plot(times, voltages, {{"label", "Node " + to_string(node_id) + "(" + ckt.node_list[node_id] + ")"}}); 
            } 
            plt::legend(); 
            plt::xlabel("Time (s)"); 
            plt::ylabel("Voltage (V)"); 
            plt::title("Transient Analysis Node Voltages"); 
            plt::grid(true); 
            plt::show();
        }
        else if (analysis.type == "HB"){
            double freq = analysis.parameters.count("freq") ? analysis.parameters["freq"] : 1e3;
            double harm = analysis.parameters.count("harm") ? analysis.parameters["harm"] : 5;
            sol.PSS_solve_harmonic_balance(analysis);
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

            // 我超级想知道为什么这段代码在solver.cpp里画不出来图
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

            std::cout << "Plotted HB time-domain voltages for all nodes.\n";
        }
        else if (analysis.type == "SHOOTING"){
            double freq;
            if (analysis.parameters.count("freq")){
                freq = analysis.parameters["freq"];
            } else {
                bool flag_found = false;
                for (auto& source : ckt.sources){
                    if (source.type == "V_SIN" || source.type == "I_SIN"){
                        if (source.parameters.count("freq")){
                            freq = source.parameters["freq"];
                            flag_found = true;
                            break;
                        }
                    }
                }
                if (!flag_found){
                    cerr << "No frequency specified for SHOOTING analysis and no sinusoidal source found.\n";
                    continue;
                }
            }
            double period_T = 1.0 / freq;
            double tstep = period_T / 1000.0; //默认时间步长为周期的1/1000
            sol.PSS_solve_shooting(period_T, tstep);
            //打印要观察的节点电压
            for (int node_id : sol.get_plot_node_ids()){
                sol.print_node_voltage(node_id);
            }
            //跑一次瞬态以获得时域波形
            // To DO: 需要一个带初值的TRAN求解接口
            sol.TRAN_solve(period_T, tstep);
            //绘制节点时域波形
            plt::figure();
            for (const auto& tran_plot_pair : sol.get_tran_plot_data()){ 
                int node_id = tran_plot_pair.first; 
                const auto& time_voltage_series = tran_plot_pair.second; 
                vector<double> times, voltages; 
                for (const auto& tv : time_voltage_series){ 
                    times.push_back(tv.first); 
                    voltages.push_back(tv.second); 
                } 
                plt::plot(times, voltages, {{"label", "Node " + to_string(node_id) + "(" + ckt.node_list[node_id] + ")"}}); 
            } 
            plt::legend(); 
            plt::xlabel("Time (s)"); 
            plt::ylabel("Voltage (V)"); 
            plt::title("Transient Analysis Node Voltages"); 
            plt::grid(true); 
            plt::show();
        }
        else {
            cerr << "Unknown analysis type: " << analysis.type << "\n";
        }

    }
    return 0;
}