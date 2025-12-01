#pragma once
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <Eigen/Dense>

class solver {
    private:
        circuit ckt;
        analysis analysis_type;
        Eigen::VectorXd J; //电流源向量
        //经过主元置换后的J向量
        Eigen::VectorXd J_permuted;
        //线性MNA矩阵
        Eigen::MatrixXd liner_Y;
        //LU分解矩阵
        Eigen::MatrixXd L;
        Eigen::MatrixXd U;
        //最终MNA矩阵
        Eigen::MatrixXd MNA_Y;
        //各节点电压
        Eigen::VectorXd node_voltages;
        //支路电流（用于电压源等）
        Eigen::VectorXd branch_currents;
        //构建线性MNA矩阵和J向量
        void build_linear_MNA();
        //非线性器件的处理
        void build_nonlinear_MNA();
        //电源的处理
        void build_sources_MNA();

        //高斯消去法线性MNA方程求解
        void solve_linear_MNA_Gauss();
        //库函数LU分解法求解MNA方程
        void solve_linear_MNA_LU();
        //手动得到LU分解矩阵
        void get_linear_MNA_LU_manual();
        //使用已计算的L和U矩阵求解方程
        void solve_with_LU_matrices();
        //Gauss-Jacobi迭代法求解线性MNA方程
        void solve_linear_MNA_Gauss_Jacobi();
        
        //MNA矩阵操作辅助函数
        void addToY(int rowNode, int colNode, double val);
        void addToJ(int node, double val);
    public:
        solver(circuit& ckt_, analysis& analysis_);
        //直流分析
        void DC_solve();
        void DC_solve(const Eigen::VectorXd& initial_voltages);
        void DC_solve(const std::map<std::string, double>& node_voltage_map);
        //瞬态分析
        void TRAN_solve();
        //稳态分析
        void Pss_solve();
};