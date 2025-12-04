#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "circuit.hpp"   // 包含 circuit/device/model 定义
#include "parse_netlist.hpp"
#include "solver.hpp"
#include <map>

using namespace std;

int main_tran(){
    circuit ckt;
    vector<string> netlist;
    vector<analysis> analyses;
    //parseNetlistFile("resistor_net.sp", ckt, analyses); //测试电阻网络
    parseNetlistFile("RLC.sp", ckt, analyses); 
    solver sol(ckt, analyses[0]);

    //设置不同的直流分析方法，测试瞬态分析
    cout << "=== 测试不同的线性方程求解方法 ===\n";
    
    // 1. 使用LU分解法（默认）
    cout << "\n1. LU分解法：\n";
    sol.setLinearSolverMethod(LinearSolverMethod::LU_DECOMPOSITION);
    sol.TRAN_solve(); // tstop=5s, tstep=0.01s
    
    // // 2. 使用高斯消去法
    // cout << "\n2. 高斯消去法：\n";
    // sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_ELIMINATION);
    // sol.TRAN_solve(5,0.01); // tstop=5s, tstep=0.01s
    
    // // 3. 使用手动LU分解法
    // cout << "\n3. 手动LU分解法：\n";
    // sol.setLinearSolverMethod(LinearSolverMethod::MANUAL_LU);
    // sol.TRAN_solve(5,0.01); // tstop=5s, tstep=0.01s
    
    // // 4. 使用Gauss-Jacobi迭代法
    // cout << "\n4. Gauss-Jacobi迭代法：\n";
    // sol.setLinearSolverMethod(LinearSolverMethod::GAUSS_JACOBI);
    // sol.TRAN_solve(5,0.01); // tstop=5s, tstep=0.01s
    

}