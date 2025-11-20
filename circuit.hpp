#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <map>


struct circuit
{
    std::unordered_map<std::string, int> node_map; //把节点名跟节点编号对应起来
    std::vector<std::string> node_list;

    int getNodeID(const std::string &name); //根据节点名获取节点编号
};
