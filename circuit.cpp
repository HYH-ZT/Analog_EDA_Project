#include "circuit.hpp"
#include <vector>
#include <string>
using namespace std;
// int circuit::getNodeID(const string &name) {
//         string nm = name;
//         if (nm=="") nm="0";
//         if (nm == "0" || nm == "gnd") { // canonical ground
//             if (node_map.find("0")==node_map.end()) {
//                 node_map["0"] = 0;
//                 if (node_list.size() == 0) node_list.push_back("0");
//             }
//             return 0;
//         }
//         auto it = node_map.find(nm);
//         if (it != node_map.end()) return it->second;
//         int id = (int)node_map.size();
//         node_map[nm] = id;
//         if (id >= (int)node_list.size()) node_list.resize(id+1);
//         node_list[id] = nm;
//         return id;
// }

// 实现 circuit::getNodeID
int circuit::getNodeID(const std::string &name) {
    if (node_map.count(name)) return node_map[name];

    int newID = node_list.size();
    node_list.push_back(name);
    node_map[name] = newID;
    return newID;
}
const model* circuit::findModelConst(const std::string& modelName) {
    for (const auto &m : models) {
        if (m.name == modelName) return &m;
    }
    return nullptr;
}

void circuit::extract_MOS_capacitances() {
    std::vector<device> new_linear_devices;
    for (const auto& dev : nonlinear_devices) {
        if (dev.type == "MOS") {
            //Cgs = Cox * W * L / 2
            //Cgd = Cox * W * L / 2
            //Cd = Cs = CJ0;
            double W = dev.parameters.at("W");
            double L = dev.parameters.at("L");

            //从模型文件中获取参数 （为什么不直接存在设备参数中？）
            const model* pmodel = findModelConst(dev.model);
            double Cox = 0.0;
            double CJ0 = 0.0;
            if (pmodel) {
                if (pmodel->parameters.count("COX")) {
                    Cox = pmodel->parameters.at("COX");
                }
                if (pmodel->parameters.count("CJ0")) {
                    CJ0 = pmodel->parameters.at("CJ0");
                }
            }
            double Cgs = Cox * W * L / 2.0;
            double Cgd = Cox * W * L / 2.0;
            double Cd = CJ0;
            double Cs = CJ0;
            // 创建线性电容器件并加入列表
            if (Cgs > 0) {
                device cgs_dev;
                cgs_dev.name = "Cgs_" + dev.name;
                cgs_dev.type = "C";
                cgs_dev.node_names = {dev.node_names[1], dev.node_names[2]}; // Gate-Source
                cgs_dev.nodes = {dev.nodes[1], dev.nodes[2]};
                cgs_dev.parameters["value"] = Cgs;
                new_linear_devices.push_back(cgs_dev);
            }
            if (Cgd > 0) {
                device cgd_dev;
                cgd_dev.name = "Cgd_" + dev.name;
                cgd_dev.type = "C";
                cgd_dev.node_names = {dev.node_names[0], dev.node_names[1]}; // Drain-Gate
                cgd_dev.nodes = {dev.nodes[0], dev.nodes[1]};
                cgd_dev.parameters["value"] = Cgd;
                new_linear_devices.push_back(cgd_dev);
            }
            if (Cd > 0) {
                device cd_dev;
                cd_dev.name = "Cd_" + dev.name;
                cd_dev.type = "C";
                cd_dev.node_names = {dev.node_names[0], "0"}; // Drain-Ground
                cd_dev.nodes = {dev.nodes[0], 0};
                cd_dev.parameters["value"] = Cd;
                new_linear_devices.push_back(cd_dev);
            }
            if (Cs > 0) {
                device cs_dev;
                cs_dev.name = "Cs_" + dev.name;
                cs_dev.type = "C";
                cs_dev.node_names = {dev.node_names[2], "0"}; // Source-Ground
                cs_dev.nodes = {dev.nodes[2], 0};
                cs_dev.parameters["value"] = Cs;
                new_linear_devices.push_back(cs_dev);
            }
        }
    }
    // 将提取的线性电容器件加入 linear_devices 列表
    linear_devices.insert(linear_devices.end(), new_linear_devices.begin(), new_linear_devices.end());
}