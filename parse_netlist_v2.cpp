// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <algorithm>
// #include <memory>
// #include "circuit_v2.hpp"

// // 去除前后空格
// static inline std::string trim(const std::string& s) {
//     size_t b = s.find_first_not_of(" \t\r\n");
//     if (b == std::string::npos) return "";
//     size_t e = s.find_last_not_of(" \t\r\n");
//     return s.substr(b, e-b+1);
// }

// bool isComment(const std::string& s) {
//     return s.empty() || s[0]=='*';
// }

// bool isModel(const std::string &line) {
//     return line.size() > 6 && line.substr(0,6) == ".MODEL";
// }

// bool isAnalysis(const std::string &line) {
//     return !line.empty() && line[0]=='.' && !isModel(line);
// }



// // 解析 .MODEL
// static void parseModelLine(const std::string& line, circuit& ckt) {
//     std::istringstream iss(line);
//     std::string dot, mname, mtype;
//     iss >> dot >> mname;
//     //iss >> mtype;

//     model md;
//     md.name = mname;
//     //md.type = mtype;

//     std::string key;
//     double val;
//     while (iss >> key >> val) {
//         md.parameters[key] = val;
//     }

//     ckt.models.push_back(md);
// }

// // 解析器件（R, C, L, V, I, M）
// static void parseDeviceLine(const std::string& line, circuit& ckt) {
//     std::istringstream iss(line);
//     //device dev;
//     //dev.rawline = line;
//     std::string name;
//     iss >> name;
//     char c = toupper(name[0]);

//     /** -------------------------
//      *  R, C, L
//      *  ------------------------- */
//     if (c=='R' || c=='C' || c=='L') {
//         std::string type = std::string(1, c);
//         std::string n1, n2;
//         double val;
//         iss >> n1 >> n2 >> val;
//         auto dev = std::make_unique<RCL>(name, type,
//                                         std::vector<std::string>{n1, n2},
//                                         std::vector<int>{ckt.getNodeID(n1), ckt.getNodeID(n2)},
//                                         val);
//         // dev.node_names = {n1, n2};
//         // dev.nodes = { ckt.getNodeID(n1), ckt.getNodeID(n2) };
//         // dev.parameters["value"] = val;
//     }

//     /** -------------------------
//      *  V, I
//      *  ------------------------- */
//     else if (c=='V' || c=='I') {
//         std::string type = std::string(1, c);

//         std::string n1, n2;
//         iss >> n1 >> n2;
//         // dev.node_names = {n1, n2};
//         // dev.nodes = { ckt.getNodeID(n1), ckt.getNodeID(n2) };

//         std::string form;
//         if (!(iss >> form)) {
//             type += "_DC"; // 默认直流
//         }
//         else if (form == "DC") {
//             double v;
//             iss >> v;
//             type += "_DC";
//             dev.parameters["DC"] = v;
//         }
//         else if (form == "SIN") {
//             double v0, amp, freq, phase;
//             iss >> v0 >> amp >> freq >> phase;
//             dev.type += "_SIN";
//             dev.parameters["DC"]   = v0;
//             dev.parameters["AMP"]  = amp;
//             dev.parameters["FREQ"] = freq;
//             dev.parameters["PHASE"]= phase;
//         }
//         else {
//             dev.type += "_UNKNOWN";
//         }
//     }


//     /** -------------------------
//      *  MOS (Mxxx ...)
//      *  格式：
//      *   Mname D G S B model W L (others)
//      *  ------------------------- */
//     else if (c=='M') {
//         dev.type = "MOS";

//         std::string D, G, S;
//         std::string modelName;
//         std::string channelType;

//         double W, L;

//         iss >> D >> G >> S >> channelType >> W >> L >> modelName;

//         dev.node_names = {D, G, S};
//         dev.nodes = {
//             ckt.getNodeID(D),
//             ckt.getNodeID(G),
//             ckt.getNodeID(S)
//         };

//         dev.model = modelName;
//         dev.parameters["W"] = W;
//         dev.parameters["L"] = L;
//         dev.parameters["TYPE"] = (channelType == "n") ? 1.0 : -1.0;

//         // 其它参数
//         std::string key;
//         double val;
//         while (iss >> key >> val) {
//             dev.parameters[key] = val;
//         }
//     }

//     else {
//         dev.type = "UNKNOWN";
//     }

//     ckt.devices.push_back(dev);
// }

// // 解析分析语句
// static void parseAnalysisLine(const std::string& line, std::vector<analysis>& analysis_list) {
//     analysis a;
//     a.type = "";
//     std::istringstream iss(line);
//     std::string cmd;
//     iss >> cmd;

//     if (cmd == ".hb") {
//         a.type = "HB";
//         double freq, harm;
//         iss >> freq >> harm;
//         a.parameters["freq"] = freq;
//         a.parameters["harm"] = harm;
//     }
//     else if (cmd == ".tran") {
//         a.type = "TRAN";
//         double step, stop;
//         iss >> step >> stop;
//         a.parameters["step"] = step;
//         a.parameters["stop"] = stop;
//     }
//     else if (cmd == ".dc") {
//         a.type = "DC";
//     }

//     analysis_list.push_back(a);
// }


// /* ==========================================================
//  *       主解析入口：parseNetlistFile()
//  * ========================================================== */
// void parseNetlistFile(const std::string& filename,
//                       circuit &ckt,
//                       std::vector<analysis>& analyses)
// {
//     std::ifstream fin(filename);
//     if (!fin.is_open()) {
//         std::cerr << "❌ 无法打开文件: " << filename << std::endl;
//         return;
//     }

//     std::string line;

//     while (std::getline(fin, line)) {
//         line = trim(line);

//         if (isComment(line)) continue;

//         if (isModel(line)) {
//             parseModelLine(line, ckt);
//             continue;
//         }

//         if (isAnalysis(line)) {
//             parseAnalysisLine(line, analyses);
//             continue;
//         }

//         // 解析器件
//         parseDeviceLine(line, ckt);
//     }
// }
