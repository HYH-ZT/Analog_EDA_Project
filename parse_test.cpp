#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "circuit.hpp"   // 包含 circuit/device/model 定义
#include "parse_netlist.hpp"
#include "solver.hpp"

using namespace std;


int main_parse_test() {
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;

    // ---- 3. 调用你的解析器 ----
    parseNetlistFile("dbmixer_dc.sp", ckt, analyses);
    // ---- 4. 输出解析结果用于检查 ----
    
    cout << "========== Node Mapping ==========\n";
    for (auto &p : ckt.node_map) {
        cout << "Node name: " << p.first << " -> ID: " << p.second << "\n";
    }

    cout << "\n========== Devices ==========\n";
    for (auto &d : ckt.linear_devices) {
        cout << "Device: " << d.name << "  Type: " << d.type << "\n";
        cout << "  Nodes: ";
        for (auto &n : d.node_names) cout << n << " ";
        cout << "\n  Node IDs: ";
        for (auto &id : d.nodes) cout << id << " ";
        cout << "\n  Parameters:\n";
        for (auto &kv : d.parameters) {
            cout << "    " << kv.first << " = " << kv.second << "\n";
        }
        if (!d.model.empty()) {
            cout << "  Model: " << d.model << "\n";
        }
        cout << "  Raw line: " << d.rawline << "\n\n";
    }
    for (auto &d : ckt.nonlinear_devices) {
        cout << "Device: " << d.name << "  Type: " << d.type << "\n";
        cout << "  Nodes: ";
        for (auto &n : d.node_names) cout << n << " ";
        cout << "\n  Node IDs: ";
        for (auto &id : d.nodes) cout << id << " ";
        cout << "\n  Parameters:\n";
        for (auto &kv : d.parameters) {
            cout << "    " << kv.first << " = " << kv.second << "\n";
        }
        if (!d.model.empty()) {
            cout << "  Model: " << d.model << "\n";
        }
        cout << "  Raw line: " << d.rawline << "\n\n";
    }
    for (auto &d : ckt.sources) {
        cout << "Source: " << d.name << "  Type: " << d.type << "\n";
        cout << "  Nodes: ";
        for (auto &n : d.node_names) cout << n << " ";
        cout << "\n  Node IDs: ";
        for (auto &id : d.nodes) cout << id << " ";
        cout << "\n  Parameters:\n";
        for (auto &kv : d.parameters) {
            cout << "    " << kv.first << " = " << kv.second << "\n";
        }
        if (!d.model.empty()) {
            cout << "  Model: " << d.model << "\n";
        }
        cout << "  Raw line: " << d.rawline << "\n\n";
    }
    cout << "\n========== Models ==========\n";
    for (auto &m : ckt.models) {
        cout << "Model: " << m.name << "\n";
        for (auto &kv : m.parameters) {
            cout << "    " << kv.first << " = " << kv.second << "\n";
        }
        cout << "\n";
    }

    cout << "========== Summary ==========\n";
    cout << "Total nodes: " << ckt.node_map.size() << "\n";
    cout << "Total devices: " << ckt.linear_devices.size() + ckt.nonlinear_devices.size() << "\n";
    cout << "Total models: " << ckt.models.size() << "\n";

    // ---- 5. 调用求解器进行分析 ----
    for (auto &analysis_type : analyses) {
        solver s(ckt, analysis_type);
        if (analysis_type.type == "DC") {
            s.DC_solve();
        }
        // 可以添加对其他分析类型的支持
    }

    return 0;
}
