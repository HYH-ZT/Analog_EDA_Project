#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <map>

struct device{
    std::string name;
    std::string type;   // R, C, L, V, I, etc.
    std::vector<std::string> node_names; //节点名称列表
    std::vector<int> nodes; //对于MOS管，nodes[0]:drain, nodes[1]:gate, nodes[2]:source, nodes[3]:bulk（我们估计用不到3）
    std::map<std::string, double> parameters;
    std::string model; //对于有模型的器件，比如MOS管
    std::string rawline; //原始输入行
};
struct model{
    std::string name;
    //std::string type; // NMOS, PMOS, BJT, etc.
    std::map<std::string, double> parameters;
};
struct analysis{
    std::string type; // DC, AC, TRAN, etc.
    std::map<std::string, double> parameters;
};
struct circuit
{
    std::unordered_map<std::string, int> node_map; //把节点名跟节点编号对应起来
    std::vector<std::string> node_list;
    std::vector<device> linear_devices;
    std::vector<device> nonlinear_devices;
    std::vector<device> sources;
    std::vector<model> models;

    int getNodeID(const std::string &name); //根据节点名获取节点编号
    const model* findModelConst(const std::string& modelName);
    //提取MOS管的寄生电容，转换为线性电容器件，加入linear_devices列表
    void extract_MOS_capacitances();
};
