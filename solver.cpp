#include "solver.hpp"
#include <iostream>
#include <string>
#include <vector>

solver::solver(circuit& ckt_, analysis& analysis_, 
               LinearSolverMethod lsm, 
               TransientMethod tm,
               SteadyStateMethod ssm)
    : ckt(ckt_), analysis_type(analysis_),
      linear_solver_method(lsm), transient_method(tm), steady_state_method(ssm) {
    liner_Y.resize(ckt.node_map.size() - 1, ckt.node_map.size() - 1); //不包含地节点
    liner_Y.setZero();
    J.resize(ckt.node_map.size() - 1);
    J.setZero();
    
    // // 设置默认求解方法
    // linear_solver_method = LinearSolverMethod::LU_DECOMPOSITION;
    // transient_method = TransientMethod::TRAPEZOIDAL;
    // steady_state_method = SteadyStateMethod::SHOOTING;
}


void solver::print_branch_current(int branch_index){
    if (branch_index >= 0 && branch_index < branch_currents.size()){
        std::cout << "I(" << ckt.sources[branch_index].name << ") " << branch_currents[branch_index] << " A\n";
    }
    else{
        std::cout << "Error: Branch current index " << branch_index << " is out of range.\n";
    }
}


void solver::set_initial_node_voltages(std::string node_name, double voltage) {
    int node_id = ckt.getNodeID(node_name);
    if (node_id != -1) {
        if (node_id == 0) {
            std::cout << "Warning: Cannot set voltage for ground node.\n";
        } else {
            node_voltages[node_id - 1] = voltage;
        }
    } else {
        std::cout << "Error: Node " << node_name << " not found in the circuit.\n";
    }
}


void solver::print_node_voltages() {
    std::cout << "Node Voltages:\n";
    for (const auto& pair : ckt.node_map) {
        const std::string& node_name = pair.first;
        int node_id = pair.second;
        if (node_id == 0) {
            std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
        } else {
            std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
        }
    }
}


void solver::print_node_voltage(const int node_id) {
    std::string node_name = ckt.node_list[node_id];
    if (node_id == 0) {
        std::cout << "Node " << node_name << " (ID " << node_id << "): 0 V (Ground)\n";
    } else if (node_id > 0 && node_id < ckt.node_map.size()) {
        std::cout << "Node " << node_name << " (ID " << node_id << "): " << node_voltages[node_id - 1] << " V\n";
    } else {
        std::cout << "Error: Node ID " << node_id << " is out of range.\n";
    }
}


void solver::print_node_voltage(const std::string& node_name) {
    int node_id = ckt.getNodeID(node_name);
    if (node_id == -1) {
        std::cout << "Error: Node " << node_name << " not found in the circuit.\n";
    } else {
        print_node_voltage(node_id);
    }
}

//解析打印变量，填充ckt.print_node_ids，以及sources中器件的printI标志

void solver::parse_print_variables() {
    ckt.print_node_ids.clear();
    ckt.plot_node_ids.clear();
    //Debug: 输出要解析的打印变量
    std::cout << analysis_type.print_variables.size() << " print variables to parse:\n";
    for (const auto& var : analysis_type.print_variables) {
        std::cout << "  " << var << "\n";
    }
    std::cout << analysis_type.plot_variables.size() << " plot variables to parse:\n";
    for (const auto& var : analysis_type.plot_variables) {
        std::cout << "  " << var << "\n";
    }

    for (const auto& var : analysis_type.print_variables) {
        if (var[0] == 'V' && var[1] == '(' && var.back() == ')') {
            //节点电压变量
            std::string node_name = var.substr(2, var.size() - 3);
            int node_id = ckt.getNodeID(node_name);
            if (node_id != -1) {
                ckt.print_node_ids.push_back(node_id);
            } else {
                std::cout << "Warning: Node " << node_name << " not found for printing voltage.\n";
            }
        } else if (var[0] == 'I' && var[1] == '(' && var.back() == ')') {
            //支路电流变量
            std::string device_name = var.substr(2, var.size() - 3);
            bool found = false;
            for (auto& src : ckt.sources) {
                if (src.name == device_name) {
                    src.printI = true; //标记该器件需要打印电流
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "Warning: Source " << device_name << " not found for printing current.\n";
            }
        } else {
            std::cout << "Warning: Unrecognized print variable format: " << var << "\n";
        }
    }
    for (const auto& var : analysis_type.plot_variables) {
        if (var[0] == 'I' && var[1] == '(' && var.back() == ')') {
            //支路电流变量
            std::string device_name = var.substr(2, var.size() - 3);
            bool found = false;
            for (auto& src : ckt.sources) {
                if (src.name == device_name) {
                    src.plotI = true; //标记该器件需要打印电流
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cout << "Warning: Source " << device_name << " not found for printing current.\n";
            }
        }
        else{
            //节点电压变量
            std::string node_name = var;
            int node_id = ckt.getNodeID(node_name);
            if (node_id != -1) {
                ckt.plot_node_ids.push_back(node_id);
            } else {
                std::cout << "Warning: Node " << node_name << " not found for printing voltage.\n";
            }
        }
    }
}

