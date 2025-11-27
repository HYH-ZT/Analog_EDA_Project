#pragma once
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <Eigen/Dense>

class solver {
    private:
        circuit ckt;
        analysis analysis_type;
        Eigen::VectorXd J; //电流源向量
        //线性MNA矩阵
        Eigen::MatrixXd liner_Y;
        //最终MNA矩阵
        Eigen::MatrixXd MNA_Y;
        //各节点电压
        std::vector<double> node_voltages;
    public:
        solver(circuit& ckt_, analysis& analysis_);
        //直流分析
        void DC_solve();
        //瞬态分析
        void TRAN_solve();
        //稳态分析
        void Pss_solve();
};