#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "circuit.hpp"   // 包含 circuit/device/model 定义
#include "parse_netlist.hpp"

using namespace std;


int main() {
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;

    // ---- 3. 调用你的解析器 ----
    parseNetlistFile("buffer.sp", ckt, analyses);
    // ---- 4. 输出解析结果用于检查 ----
    
    cout << "========== Node Mapping ==========\n";
    for (auto &p : ckt.node_map) {
        cout << "Node name: " << p.first << " -> ID: " << p.second << "\n";
    }

    cout << "\n========== Devices ==========\n";
    for (auto &d : ckt.liner_devices) {
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
    for (auto &d : ckt.nonliner_devices) {
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
    cout << "Total devices: " << ckt.liner_devices.size() + ckt.nonliner_devices.size() << "\n";
    cout << "Total models: " << ckt.models.size() << "\n";

    return 0;
}
