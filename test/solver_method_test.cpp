#include <iostream>
#include "circuit.hpp"
#include "parse_netlist.hpp" 
#include "solver.hpp"

using namespace std;

int main_dc2(){
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;
    
    // 解析电路网表
    parseNetlistFile("CSamp.sp", ckt, analyses); 
    
    // 创建求解器
    solver sol(ckt, analyses[0]);
    
    // 测试不同的线性方程求解方法
    cout << "=== 测试不同的线性方程求解方法 ===\n";
    
    // 1. 使用LU分解法（默认）
    cout << "\n1. LU分解法：\n";
    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
    sol.DC_solve();
    
    // 2. 使用高斯消去法
    cout << "\n2. 高斯消去法：\n";
    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_ELIMINATION);
    sol.DC_solve();
    
    // 3. 使用手动LU分解法
    cout << "\n3. 手动LU分解法：\n";
    sol.setLinearSolverMethod(LinearSolverMethod::MANUAL_LU);
    sol.DC_solve();
    
    // 4. 使用Gauss-Jacobi迭代法
    cout << "\n4. Gauss-Jacobi迭代法：\n";
    sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_JACOBI);
    sol.DC_solve();
    
    // 设置瞬态分析方法
    cout << "\n=== 设置瞬态分析方法 ===\n";
    sol.setTransientMethod(TransientMethod::TRAPEZOIDAL);
    cout << "瞬态分析方法设置为：梯形积分法\n";
    
    // 设置稳态分析方法
    cout << "\n=== 设置稳态分析方法 ===\n"; 
    sol.setSteadyStateMethod(SteadyStateMethod::SHOOTING);
    cout << "稳态分析方法设置为：Shooting Method\n";
    
    // 显示当前设置
    cout << "\n=== 当前求解方法设置 ===\n";
    cout << "线性求解方法: " << (int)sol.getLinearSolverMethod() << "\n";
    cout << "瞬态分析方法: " << (int)sol.getTransientMethod() << "\n"; 
    cout << "稳态分析方法: " << (int)sol.getSteadyStateMethod() << "\n";
    
    return 0;
}