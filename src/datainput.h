/********************************************************************************************************
* datainput.hh：代码参数解析器的实现
*********************************************************************************************************/

#ifndef DATAINPUT_H_
#define DATAINPUT_H_

#include <iostream>
#include <iterator>
#include <list>
#include <string>
#include <vector>

#include "classes.h"
#include "datainput.h"

namespace MITHRA
{

  /*这类函数用于读取参数的文本文件并将其解析为大流士
  *解决者。*/
  class ParseDarius
  {

  private:

    /*解析值所需的参数和数据文件。*/
    std::list<std::string>&		jobFile_;
    Mesh& 				mesh_;
    Bunch&				bunch_;
    Seed&				seed_;
    std::vector<Undulator>&		undulator_;
    std::vector<ExtField>&              extField_;
    std::vector<FreeElectronLaser>&	FEL_;

  public:

    ParseDarius 	(std::list<std::string>& jobFile, Mesh& mesh, Bunch& bunch, Seed& seed,
			 std::vector<Undulator>& undulator, std::vector<ExtField>& extField,
			 std::vector<FreeElectronLaser>& FEL);

    /*从文件中读取参数，并为FEL仿真设置所有解析参数。*/
    void setJobParameters ();

    /*在求解器中读取为网格解析的参数。*/
    void readMesh 	(std::list <std::string>::iterator & iter);

    /*读取darius解算器中为一组解析的参数。*/
    void readBunch 	(std::list <std::string>::iterator & iter);

    /*读取darius解算器中为种子解析的参数。*/
    void readField 	(std::list <std::string>::iterator & iter);

    /*在求解器中读取为网格解析的参数。*/
    void readUndulator 	(std::list <std::string>::iterator & iter);

    /*读取darius解算器中为种子解析的参数。*/
    void readExtField 	(std::list <std::string>::iterator & iter);

    /*读取在darius解算器中为FEL输出解析的参数。*/
    void readFEL 	(std::list <std::string>::iterator & iter);
  };
}
#endif
