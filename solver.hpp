#pragma once
#include "circuit.hpp"
#include "parse_netlist.hpp"
#include <Eigen/Dense>

// 线性方程求解方法枚举
enum class LinearSolverMethod {
    GAUSS_ELIMINATION = 0,  // 高斯消去法
    LU_DECOMPOSITION = 1,   // LU分解法（库函数）
    MANUAL_LU = 2,          // 手动LU分解法
    GAUSS_JACOBI = 3        // Gauss-Jacobi迭代法
};

// 瞬态分析方法枚举
enum class TransientMethod {
    FORWARD_EULER = 0,      // 前向欧拉法
    BACKWARD_EULER = 1,     // 后向欧拉法（隐式）
    TRAPEZOIDAL = 2         // 梯形积分法
};

// 稳态分析方法枚举
enum class SteadyStateMethod {
    NEWTON_RAPHSON = 0,     // Newton-Raphson法
    DAMPED_NEWTON = 1,      // 阻尼Newton法
    CONTINUATION = 2        // 延拓法
};

class solver {
    private:
        circuit ckt;
        analysis analysis_type;
        
        // 求解方法选择
        LinearSolverMethod linear_solver_method;
        TransientMethod transient_method;
        SteadyStateMethod steady_state_method;
        
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
        //动态器件名称到支路电流索引的映射
        std::map<std::string, int> dynamic_device_current_map;
        //构建线性MNA矩阵和J向量
        void build_linear_MNA(bool in_tran = false);
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
        //根据设置的方法选择线性方程求解算法
        void solve_linear_MNA();
        //根据输入参数，选择矩阵方程求解方法（兼容旧接口）
        //method: 0-高斯消去法，1-LU分解法，2-手动LU分解法，3-Gauss-Jacobi迭代法
        void solve_linear_MNA(int method);
        
        //MNA矩阵操作辅助函数
        void addToY(int rowNode, int colNode, double val);
        void addToJ(int node, double val);
        //解析打印变量，填充ckt.print_node_ids
        void parse_print_variables();

        //构建瞬态分析电路
        void build_transient_ckt(double tstep);


    public:
        solver(circuit& ckt_, analysis& analysis_);
        
        // 设置求解方法
        void setLinearSolverMethod(LinearSolverMethod method) { linear_solver_method = method; }
        void setTransientMethod(TransientMethod method) { transient_method = method; }
        void setSteadyStateMethod(SteadyStateMethod method) { steady_state_method = method; }
        
        // 获取当前求解方法
        LinearSolverMethod getLinearSolverMethod() const { return linear_solver_method; }
        TransientMethod getTransientMethod() const { return transient_method; }
        SteadyStateMethod getSteadyStateMethod() const { return steady_state_method; }
        
        //直流分析
        void DC_solve();
        void DC_solve(const Eigen::VectorXd& initial_voltages,bool in_tran = false);
        void DC_solve(const std::map<std::string, double>& node_voltage_map,bool in_tran = false);
        //瞬态分析
        void TRAN_solve();
        //
        void TRAN_solve(double tstop, double tstep);
        //稳态分析

};