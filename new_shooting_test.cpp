#include <iostream>
#include <vector>
#include <string>
#include <matplotlibcpp.h>
#include <Python.h>
#include "solver.hpp"
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <cmath>

#ifndef M_PI
constexpr double M_PI = 3.141592653589793238462643383279502884;
#endif

using namespace std;
namespace plt = matplotlibcpp;


#include <unordered_map>


int main(){
    string netlist_file;
    cin >> netlist_file;
    circuit ckt;
    vector<analysis> analyses;
    parseNetlistFile(netlist_file, ckt, analyses);
    if (analyses.empty()) {
        cerr << "No analysis specified in the netlist.\n";
        return 1;
    }

    solver sol(ckt, analyses[0]);
    double freq;
    if (analyses[0].parameters.count("freq")){
        freq = analyses[0].parameters["freq"];
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
        }
    }

    //输入采样点数
    cout << "Enter number of samples for SHOOTING analysis (e.g., 1000): ";
    double N_sample;
    cin >> N_sample;

    double period_T = 1.0 / freq;
    double tstep = period_T / N_sample; //默认时间步长为周期的1/1000
    sol.PSS_solve_shooting_new_new(period_T, tstep, 100, 1e-9);

    // double tstep = analyses[0].parameters.at("tstep");
    // double tstop = analyses[0].parameters.at("tstop");
    // sol.TRAN_solve_new_new(tstop, tstep);

    // cout << "\nbegin to draw figure.\n\n";

    // plot_tran_results(ckt, sol);

    // for (const auto& tran_plot_pair : sol.get_tran_plot_data()){
    //     int node_id = tran_plot_pair.first; 
    //     const auto& time_voltage_series = tran_plot_pair.second; 
    //     for (const auto& tv : time_voltage_series){ 
    //         cout << tv.first << "\t" << tv.second << "\n"; 
    //     } 
    //     system("pause");
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
    return 0;
}
