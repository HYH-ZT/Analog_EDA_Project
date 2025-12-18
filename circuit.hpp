#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <Eigen/Dense>

struct device{
    std::string name;
    std::string type;   // R, C, L, V, I, etc.
    std::vector<std::string> node_names; //节点名称列表
    std::vector<int> nodes; //节点编号列表
    std::map<std::string, double> parameters;   //器件参数列表
    std::string model; //对于有模型的器件，比如MOS管
    std::string rawline; //原始输入行，调试用

    //对于电压源以及电感等需要引入支路电流变量的器件，记录其支路电流变量在MNA矩阵中的索引
    int branch_current_index = -1; 
    //如果是动态等效器件，记录原始器件名称，否则为空
    std::string original_device_name = ""; 
    //是否在分析中打印该器件电流
    bool printI = false; 

    bool plotI = false;  //是否在分析中绘制该器件电流
    //电感等效来的电压源，指向原始电感器件指针
    device* original_device_ptr = nullptr;
};
struct model{
    std::string name;
    //std::string type; // NMOS, PMOS, BJT, etc.
    std::map<std::string, double> parameters;
};
struct analysis{
    std::string type; // DC, AC, TRAN, etc.
    std::map<std::string, double> parameters;
    std::vector<std::string> print_variables; // .print命令指定的输出变量,如 V(103), I(VTN)
    std::vector<std::string> plot_variables; // .plotnv命令指定的绘图变量
};
struct circuit
{
    std::unordered_map<std::string, int> node_map; //把节点名跟节点编号对应起来
    std::vector<std::string> node_list;
    int getNodeID(const std::string &name); //根据节点名获取节点编号

    std::vector<device> linear_devices;
    std::vector<device> nonlinear_devices;
    std::vector<device> sources;
    std::vector<model> models;

    //添加PLOTNV指令需要展示的节点
    std::vector<std::string> plot_node_names;
    std::vector<int> plot_node_ids;
    std::vector<int> plot_branch_current_indices;
    //添加需要打印的node电压和sources的电流
    std::vector<int> print_node_ids;
    std::vector<int> print_branch_current_indices;

    const model* findModelConst(const std::string& modelName);
    //提取MOS管的寄生电容，转换为线性电容器件，加入linear_devices列表
    void extract_MOS_capacitances();

    //===========================================
    //shooting method专用变量与函数
    //===========================================
    int allocate_internal_node();
    void add_resistor(int n1, int n2, double value);
    int add_voltage_source(int n1, int n2, double value);
    int add_current_source(int n1, int n2, double value);
};


//==========================================
// 专用于shooting method的电路结构
//==========================================

//电容状态结构体
struct CapacitorState {
    double v_prev = 0.0;   // 上一步电容电压
    double i_prev = 0.0;   // 上一步电容电流（可选，但建议保留）
};

//电容骨架结构体
struct CapacitorSkeleton {
    int n1, n2;        // 连接的两个节点
    int mid;           // 内部节点编号
    double Req;        // 固定等效电阻
    int veq_index;     // 在 sources 中的 index
    int state_index = -1;   // ★ 必须
};

//电感状态结构体
struct InductorState {
    double i_prev;     // 上一步电感电流
    double v_prev;     // 上一步电感电压（梯形法需要）
};

//电感骨架结构体
struct InductorSkeleton {
    int n1, n2;        // 原始端点
    double L;
    double Req;        // 2L / dt
    int ieq_index;     // 等效电流源索引
    int state_index = -1;   // ★ 必须
};