/********************************************************************************************************
* readdata.hh：读取作业文件行函数的实现。
********************************************************************************************************/

#ifndef READDATA_HH_
#define READDATA_HH_

#include <list>
#include <string>
#include <vector>

#include "datainput.h"
#include "fieldvector.h"
#include "stdinclude.h"

namespace MITHRA
{

  /*从输入文件中读取数据并将其存储到字符串列表中。*/
  std::list <std::string> 	read_file 		(char const * filename);

  /*清理存储的字符串向量并使其有组织。*/
  void 				cleanJobFile 		(std::list <std::string> & jobFile);

  /*读取该行的ParamaterName。*/
  std::string 			parameterName 		(std::string line);

  /*读取字符串参数的值。*/
  std::string 			stringValue 		(std::string line);

  /*双参数的读值。*/
  Double 			doubleValue 		(std::string line);

  /*读取整型参数的值。*/
  int 				intValue 		(std::string line);

  /*读取布尔参数的值。*/
  bool 				boolValue 		(std::string line);

  /*读取一个矢量参数的值。*/
  std::vector <Double> 		vectorDoubleValue 	(std::string line);

  /*读取一个矢量参数的值。*/
  std::vector <unsigned int> 	vectorIntValue 		(std::string line);

  /*map参数的读值。*/
  void 				mapValue 		(std::string line, unsigned int & tag, std::string & model);
}
#endif
