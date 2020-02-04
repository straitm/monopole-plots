#include "Data_Quality.hh"

#include "Constants.hh"

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>



Data_Quality::Data_Quality(std::string line)
{
  std::vector<std::string> cells;
  std::stringstream line_stream(line);
  std::string cell;

  while (std::getline(line_stream, cell, ','))
    cells.push_back(cell);

  // std::cout << "number of identified cells: " << cells.size()
  // 	    << std::endl;

  // unsigned n_cell = 0;
  // for (auto const& c : cells)
  //   std::cout << n_cell++ << " " << c << std::endl;

  run_             = std::stoi(cells.at(1));
  subrun_          = std::stoi(cells.at(2));
  flag_            = std::stoi(cells.at(4));
  comment_         = cells.at(5);
  n_good_diblocks_ = std::stoi(cells.at(7));
}



std::string Data_Quality::comment() const
{
  return comment_;
}



int Data_Quality::flag() const
{
  return flag_;
}



bool Data_Quality::is_good() const
{
  if (flag_ != 0)
    return false;

  return true;
}



bool Data_Quality::is_good_with_14_diblocks() const
{
  if (!is_good())
    return false;

  if (n_good_diblocks() < 14)
    return false;

  return true;
}



int Data_Quality::n_good_diblocks() const
{
  return n_good_diblocks_;
}



int Data_Quality::run() const
{
  return run_;
}



int Data_Quality::subrun() const
{
  return subrun_;
}



dq_t parse_data_quality_csv_file()
{
  dq_t result;

  std::string file_name = DATA_QUALITY_FILE;
  std::ifstream file(file_name.c_str());
  if (!file.is_open())
  {
    std::cerr << file_name << " does not exist!" << std::endl;
    assert(false);
  }

  std::string line;
  unsigned n = 0;
  while (getline(file, line))
  {
    ++n;

    if (n == 1)
    {
      std::cout << "DQParser: header = " << line << "/n/n" << std::endl;
    } else {
      // std::cout << "parsing: " << line << std::endl;
      Data_Quality info(line);

      // let's assume that the database entries are unique:
      result[info.run()][info.subrun()] = info;
   }
  }

  std::cout << "DQ Parser: " << n - 1 << " lines read in" << std::endl;
  
  return result;
}
