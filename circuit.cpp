#include "circuit.hpp"
#include <string>
using namespace std;
int circuit::getNodeID(const string &name) {
        string nm = name;
        if (nm=="") nm="0";
        if (nm == "0" || nm == "gnd") { // canonical ground
            if (node_map.find("0")==node_map.end()) {
                node_map["0"] = 0;
                if (node_list.size() == 0) node_list.push_back("0");
            }
            return 0;
        }
        auto it = node_map.find(nm);
        if (it != node_map.end()) return it->second;
        int id = (int)node_map.size();
        node_map[nm] = id;
        if (id >= (int)node_list.size()) node_list.resize(id+1);
        node_list[id] = nm;
        return id;
}