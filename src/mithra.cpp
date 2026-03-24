/********************************************************************************************************
 MITHRA (Completely Numerical Calculation of Free Electron Laser Radiation)
 Version 2.0, copyright 2019, Arya Fallahi
 Code developers:
 - Arya Fallahi ( IT'IS Foundation, Zurich, Switzerland )
 - Arnau Albe   ( Paul Scherrer Institut (PSI), Villigen, Switzerland )
 ********************************************************************************************************
 mithra.cpp: Main program file
 MPI version should be compiled in the mithra folder using the make file:
 make
 to commit to git use
 git add -u - git commit -m "message" - git push
 ********************************************************************************************************/


#include <list>
#include <mpi.h>
#include <string>
#include <sys/time.h>
#include <vector>
#include <iostream>

#include "stdinclude.h"
#include "classes.h"
#include "database.h"
#include "datainput.h"
#include "fdtd.h"
#include "fdtdSC.h"
#include "fieldvector.h"
#include "readdata.h"
#include "solver.h"

int main (int argc, char * (argv) [])
{

  /*初始化MPI，退出时自动完成finalize*/
  MPI_Init(&argc,&argv);

  /*激活名称空间*/
  using namespace MITHRA;

  /*当我们开始模拟时检索时间。*/
  timeval simulationStart, simulationEnd;
  gettimeofday(&simulationStart, NULL);

  helloMessage();

  /*解析命令行选项*/
  std::list<std::string> jobFile = read_file(argv[1]);
  cleanJobFile(jobFile);

  /*创建求解器数据库。*/
  Mesh                                  mesh;
  mesh.initialize();

  /*创建集群数据库。*/
  Bunch                                 bunch;

  /*创建种子数据库。*/
  Seed                                  seed;

  /*创建波动器数据库。*/
  std::vector<Undulator>                undulator;
  undulator.clear();

  /*创建外部字段数据库。*/
  std::vector<ExtField>                 extField;
  extField.clear();

  /*创建自由电子激光数据库。*/
  std::vector<FreeElectronLaser>        FEL;
  FEL.clear();

  /*打开输入参数解析器并实例化数据库。*/
  ParseDarius parser (jobFile, mesh, bunch, seed, undulator, extField, FEL);
  parser.setJobParameters();

  /*显示模拟的参数。*/
  mesh.show();
  bunch.show();
  seed.show();

  for (unsigned int i = 0; i < undulator.size(); i++) 	undulator[i].show();
  for (unsigned int i = 0; i < extField.size();  i++) 	extField[i] .show();

  /*初始化用于FDTD计算的类。*/
  Solver *solver;
  if ( mesh.spaceCharge_ )
    solver = new FdTdSC (mesh, bunch, seed, undulator, extField, FEL);
  else
    solver = new FdTd (mesh, bunch, seed, undulator, extField, FEL);

  /*求解指定时间内的场和束分布。*/
  solver->solve();

  /*计算总模拟时间。*/
  gettimeofday(&simulationEnd, NULL);
  Double deltaTime = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
  deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );
  printmessage(std::string(__FILE__), __LINE__, std::string("::: total simulation time [seconds] = ") + stringify(deltaTime) );

  MPI_Finalize();
}
