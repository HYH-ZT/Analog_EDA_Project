#include <iostream>
#include <vector>
#include <string>
#include "solver.hpp"
#include "circuit.hpp"
#include "parse_netlist.hpp"

using namespace std;

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
        // else if (analysis.type == "TRAN"){
        //     double tstep = analysis.parameters.count("step") ? analysis.parameters["step"] : 1e-3;
        //     double tstop = analysis.parameters.count("stop") ? analysis.parameters["stop"] : 1.0;
        //     sol.TRAN_solve(tstop, tstep);
        // }
        // else if (analysis.type == "HB"){
        //     double freq = analysis.parameters.count("freq") ? analysis.parameters["freq"] : 1e3;
        //     double harm = analysis.parameters.count("harm") ? analysis.parameters["harm"] : 5;
        //     sol.HB_solve(freq, (int)harm);
        // }
        // else {
        //     cerr << "Unknown analysis type: " << analysis.type << "\n";
        // }

    }
    return 0;
}