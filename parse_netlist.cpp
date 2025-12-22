#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "circuit.hpp"

// 去除前后空格
static inline std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e-b+1);
}

bool isComment(const std::string& s) {
    return s.empty() || s[0]=='*';
}

bool isModel(const std::string &line) {
    return line.size() > 6 && line.substr(0,6) == ".MODEL";
}

bool isAnalysis(const std::string &line) {
    return !line.empty() && line[0]=='.' && !isModel(line);
}



// 解析 .MODEL
static void parseModelLine(const std::string& line, circuit& ckt) {
    std::istringstream iss(line);
    std::string dot, mname, mtype;
    iss >> dot >> mname;
    //iss >> mtype;

    model md;
    md.name = mname;
    //md.type = mtype;

    std::string key;
    double val;
    while (iss >> key >> val) {
        md.parameters[key] = val;
    }

    ckt.models.push_back(md);
}

// 解析器件（R, C, L, V, I, M）
static void parseDeviceLine(const std::string& line, circuit& ckt) {
    std::istringstream iss(line);
    device dev;
    dev.rawline = line;

    iss >> dev.name;
    char c = toupper(dev.name[0]);

    /** -------------------------
     *  R, C, L
     *  ------------------------- */
    if (c=='R' || c=='C' || c=='L') {
        dev.type = std::string(1, c);
        std::string n1, n2;
        double val;
        iss >> n1 >> n2 >> val;

        dev.node_names = {n1, n2};
        dev.nodes = { ckt.getNodeID(n1), ckt.getNodeID(n2) };
        dev.parameters["value"] = val;
    }

    /** -------------------------
     *  V, I
     *  ------------------------- */
    else if (c=='V' || c=='I') {
        dev.type = std::string(1, c);

        std::string n1, n2;
        iss >> n1 >> n2;
        dev.node_names = {n1, n2};
        dev.nodes = { ckt.getNodeID(n1), ckt.getNodeID(n2) };

        std::string form;
        if (!(iss >> form)) {
            dev.type += "_DC"; // 默认直流
        }
        else if (form == "DC") {
            double v;
            iss >> v;
            dev.type += "_DC";
            dev.parameters["DC"] = v;
        }
        else if (form == "SIN") {
            double v0, amp, freq, phase;
            iss >> v0 >> amp >> freq >> phase;
            dev.type += "_SIN";
            dev.parameters["DC"]   = v0;
            dev.parameters["AMP"]  = amp;
            dev.parameters["FREQ"] = freq;
            dev.parameters["PHASE"]= phase;
        }
        else {
            dev.type += "_UNKNOWN";
        }
    }


    /** -------------------------
     *  MOS (Mxxx ...)
     *  格式：
     *   Mname D G S B model W L (others)
     *  ------------------------- */
    else if (c=='M') {
        dev.type = "MOS";

        std::string T1, G, T2;
        std::string modelName;
        std::string channelType;

        double W, L;

        iss >> T1 >> G >> T2 >> channelType >> W >> L >> modelName;

        dev.node_names = {T1, G, T2};
        dev.nodes = {
            ckt.getNodeID(T1),
            ckt.getNodeID(G),
            ckt.getNodeID(T2)
        };

        dev.model = modelName;
        dev.parameters["W"] = W;
        dev.parameters["L"] = L;
        dev.parameters["TYPE"] = (channelType == "n") ? 1.0 : -1.0;

        // 其它参数
        std::string key;
        double val;
        while (iss >> key >> val) {
            dev.parameters[key] = val;
        }
    }

    else {
        dev.type = "UNKNOWN";
    }

    if (dev.type == "MOS") {
        ckt.nonlinear_devices.push_back(dev);
    } else if (c == 'V' || c == 'I') {
        ckt.sources.push_back(dev);
    } else {
        ckt.linear_devices.push_back(dev);
    }
    //ckt.devices.push_back(dev);
}

// 解析分析语句
static void parseAnalysisLine(const std::string& line, std::vector<analysis>& analysis_list, circuit& ckt) {
    analysis a;
    a.type = "";
    std::istringstream iss(line);
    std::string cmd;
    iss >> cmd;

    if (cmd == ".hb" || cmd == ".HB") {
        a.type = "HB";
        double freq, harm;
        iss >> freq >> harm;
        a.parameters["freq"] = freq;
        a.parameters["harm"] = harm;
        analysis_list.push_back(a);
    }
    else if (cmd == ".tran" || cmd == ".TRAN") {
        a.type = "TRAN";
        double step, stop;
        iss >> step >> stop;
        a.parameters["tstep"] = step;
        a.parameters["tstop"] = stop;
        analysis_list.push_back(a);
    }
    else if (cmd == ".dc") {
        a.type = "DC";
        analysis_list.push_back(a);
    }
    else if (cmd == ".shooting" || cmd == ".SHOOTING") {
        a.type = "SHOOTING";
        double freq;
        if (iss >> freq) {
            a.parameters["freq"] = freq;
        }
        // else：什么都不做，表示 freq 未指定
        analysis_list.push_back(a);
    }
    else if (cmd == ".print") {
        // .print命令: .print tran V(103) V(105) I(VTN) I(VTP)
        std::string analysis_type;
        iss >> analysis_type; // 读取分析类型(tran, dc, ac等)
        
        // 将分析类型转换为大写
        for (auto& c : analysis_type) {
            c = std::toupper(static_cast<unsigned char>(c));
        }
        
        // 查找最近的匹配分析类型
        bool found = false;
        for (int i = analysis_list.size() - 1; i >= 0; --i) {
            if (analysis_list[i].type == analysis_type) {
                // 读取所有输出变量并添加到该分析
                std::string var;
                while (iss >> var) {
                    analysis_list[i].print_variables.push_back(var);
                }
                found = true;
                break;
            }
        }
        
        if (!found) {
            std::cerr << "⚠️  警告: .print命令未找到对应的 " << analysis_type << " 分析" << std::endl;
        }
    }
    else if (cmd == ".plotnv" || cmd == ".PLOTNV") {
        // .PLOTNV命令: .PLOTNV V(1) V(2) V(3)
        std::string nodeName;
        // while (iss >> nodeName) {
        //     ckt.plot_node_names.push_back(nodeName);
        //     int nodeID = ckt.getNodeID(nodeName);
        //     ckt.plot_node_ids.push_back(nodeID);
        // }
        
        //处理电流
        std::vector<std::string> parts;
        while(iss >> nodeName){
            parts.push_back(nodeName);
        }
        //加入到所有分析中
        for(auto& a : analysis_list){
            for(const auto& part : parts){
                a.plot_variables.push_back(part);
            }
        }
    }
}


/* ==========================================================
 *       主解析入口：parseNetlistFile()
 * ========================================================== */
void parseNetlistFile(const std::string& filename,
                      circuit &ckt,
                      std::vector<analysis>& analyses)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "❌ 无法打开文件: " << filename << std::endl;
        return;
    }

    std::string line;

    ckt.getNodeID("0"); //确保节点0存在

    while (std::getline(fin, line)) {
        line = trim(line);

        if (isComment(line)) continue;

        if (isModel(line)) {
            parseModelLine(line, ckt);
            continue;
        }

        if (isAnalysis(line)) {
            parseAnalysisLine(line, analyses, ckt);
            continue;
        }

        // 解析器件
        parseDeviceLine(line, ckt);
    }
}
