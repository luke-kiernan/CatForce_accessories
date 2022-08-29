
#include "LifeAPI.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <vector>

std::string GetRLE(const std::vector<std::vector<bool>> &life2d) {
  if (life2d.empty())
    return "";

  if (life2d[0].empty())
    return "";

  std::stringstream result;

  unsigned eol_count = 0;

  for (unsigned j = 0; j < life2d[0].size(); j++) {
    bool last_val = life2d[0][j];
    unsigned run_count = 0;

    for (const auto & i : life2d) {
      bool val = i[j];

      // Flush linefeeds if we find a live cell
      if (val && eol_count > 0) {
        if (eol_count > 1)
          result << eol_count;

        result << "$";

        eol_count = 0;
      }

      // Flush current run if val changes
      if (val == !last_val) {
        if (run_count > 1)
          result << run_count;

        if (last_val == 1)
          result << "o";
        else
          result << "b";

        run_count = 0;
      }

      run_count++;
      last_val = val;
    }

    // Flush run of live cells at end of line
    if (last_val) {
      if (run_count > 1)
        result << run_count;

      result << "o";

      run_count = 0;
    }

    eol_count++;
  }

  // Flush trailing linefeeds
  if (eol_count > 0) {
    if (eol_count > 1)
      result << eol_count;

    result << "$";

    eol_count = 0;
  }

  return result.str();
}

std::string GetRLE(const LifeState &s) {
  // xybounds here is with (0,0) being at the center of the torus.
  std::array<int, 4> xyBounds = s.XYBounds();
  size_t startHoriz = size_t((xyBounds[0]+64) % 64);
  size_t startVert = size_t((xyBounds[1]+64) % 64);
  size_t width = xyBounds[2]-xyBounds[0]+1;
  size_t height = xyBounds[3]-xyBounds[1]+1;
  //std::cout << xyBounds[0] << " " << xyBounds[1] << " ";
  //std::cout << xyBounds[2] << " " << xyBounds[3] << std::endl;
  //std::cout << startHoriz << " " << startVert << std::endl;

  std::vector<std::vector<bool>> vec(width, std::vector<bool>(height));
  // a list of columns.
  /*for(unsigned j = 0; j < 63; ++j){
    for( unsigned i = 0; i < 63; ++i){
      if(s.GetCell(i,j) == 1){
        std::cout << i << " " << j <<std::endl;
        assert(startHoriz <= i && i < startHoriz + width);
        assert(startVert <= j && j < startVert + height);
      }
    }
  }*/ 

  for (unsigned j = 0; j < height; j++)
    for (unsigned  i = 0; i < width; i++)
      vec[i][j] = s.GetCell(i+startHoriz, j+startVert) == 1;
  
  return GetRLE(vec);
}

std::vector<SymmetryTransform> SymmetryChainFromEnum(const StaticSymmetry sym) {
  switch (sym) {
  case StaticSymmetry::C1:
    return {};
  case StaticSymmetry::D2AcrossY:
    return {ReflectAcrossY};
  case StaticSymmetry::D2AcrossYEven:
    return {ReflectAcrossYEven};
  case StaticSymmetry::D2AcrossX:
    return {ReflectAcrossX};
  case StaticSymmetry::D2AcrossXEven:
    return {ReflectAcrossXEven};
  case StaticSymmetry::D2diagodd:
    return {ReflectAcrossYeqX};
  case StaticSymmetry::D2negdiagodd:
    return {ReflectAcrossYeqNegXP1};
  case StaticSymmetry::C2:
    return {Rotate180OddBoth};
  case StaticSymmetry::C2even:
    return {Rotate180EvenBoth};
  case StaticSymmetry::C2horizontaleven:
    return {Rotate180EvenHorizontal};
  case StaticSymmetry::C2verticaleven:
    return {Rotate180EvenVertical};
  case StaticSymmetry::C4:
    return {Rotate90, Rotate180OddBoth};
  case StaticSymmetry::C4even:
    return {Rotate90Even, Rotate180EvenBoth};
  case StaticSymmetry::D4:
    return {ReflectAcrossX, ReflectAcrossY};
  case StaticSymmetry::D4even:
    return {ReflectAcrossXEven, ReflectAcrossYEven};
  case StaticSymmetry::D4horizontaleven:
    return {ReflectAcrossYEven, ReflectAcrossX};
  case StaticSymmetry::D4verticaleven:
    return {ReflectAcrossXEven, ReflectAcrossY};
  case StaticSymmetry::D4diag:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegXP1};
  case StaticSymmetry::D4diageven:
    return {ReflectAcrossYeqX, ReflectAcrossYeqNegX};
  case StaticSymmetry::D8:
    return {Rotate90, Rotate180OddBoth, ReflectAcrossYeqX};
  case StaticSymmetry::D8even:
    return {Rotate90Even, Rotate180EvenBoth, ReflectAcrossYeqX};
  }
}

StaticSymmetry SymmetryFromString(const std::string &name) {
  std::string start = name.substr(0, 2);
  std::string rest = name.substr(2);
  if (start == "D2") {
    if (rest == "-" or rest == "vertical") {
      return StaticSymmetry::D2AcrossX;
    } else if (rest == "-even" or rest == "verticaleven") {
      return StaticSymmetry::D2AcrossXEven;
    } else if (rest == "|" or rest == "horizontal") {
      return StaticSymmetry::D2AcrossY;
    } else if (rest == "|even" or rest == "horizontaleven") {
      return StaticSymmetry::D2AcrossYEven;
    } else if (rest == "/" or rest == "/odd") {
      return StaticSymmetry::D2negdiagodd;
    } else if (rest == "\\" or rest == "\\odd") {
      return StaticSymmetry::D2diagodd;
    }
  } else if (start == "C2") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C2;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C2even;
    } else if (rest == "horizontaleven" or rest == "|even") {
      return StaticSymmetry::C2horizontaleven;
    } else if (rest == "verticaleven" or rest == "-even" or rest == "_2") {
      return StaticSymmetry::C2verticaleven;
    }
  } else if (start == "C4") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::C4;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::C4even;
    }
  } else if (start == "D4") {
    std::string evenOddInfo = rest.substr(1);
    if (rest[0] == '+' or (rest.size() > 1 and rest[1] == '+')) {
      if (evenOddInfo == "" or rest == "_+1") {
        return StaticSymmetry::D4;
      } else if (evenOddInfo == "even" or rest == "_+4") {
        return StaticSymmetry::D4even;
      } else if (evenOddInfo == "verticaleven" or evenOddInfo == "-even" or
                 rest == "_+2") {
        return StaticSymmetry::D4verticaleven;
      } else if (evenOddInfo == "horizontaleven" or evenOddInfo == "|even") {
        return StaticSymmetry::D4horizontaleven;
      }
    } else if (rest[0] == 'x' or (rest.size() > 1 and rest[1] == 'x')) {
      if (evenOddInfo == "" or rest == "_x1") {
        return StaticSymmetry::D4diag;
      } else if (evenOddInfo == "even" or rest == "_x4") {
        return StaticSymmetry::D4diageven;
      }
    }
  } else if (start == "D8") {
    if (rest == "" or rest == "_1") {
      return StaticSymmetry::D8;
    } else if (rest == "even" or rest == "_4") {
      return StaticSymmetry::D8even;
    }
  }
  return StaticSymmetry::C1;
}

std::string FilterForRegion(LifeState & state){
  // returns: rle x y.
  std::array<int,4> xybounds = state.XYBounds();
  std::string rle = GetRLE(state);
  return rle + " " + std::to_string(xybounds[0]) + " " + std::to_string(xybounds[1]);
}

int main(int argc, char *argv[]) {
  // input: activeRegionRLE, x0, y0, symmetry, earliestGen, latestGen
  // output: proper filter input for flippers, 1 orfilter per line.


  if (argc != 7){
    std::cout << "too many or too few arguments" << std::endl;
    for ( size_t i = 0; i < argc; ++i ){
      std::cout << std::string(argv[i]) << std::endl;
    }
    std::cout << "proper input: activeRegionRLE x0 y0 symmetry firstFilterGen lastFilterGen" << std::endl;
    return 0;
  }

  int x0 = std::atoi(argv[2]);
  int y0 = std::atoi(argv[3]);

  
  LifeState activeRegion = LifeState::Parse(argv[1], x0, y0, Identity);
  StaticSymmetry sym = SymmetryFromString(std::string(argv[4]));
  auto symChain = SymmetryChainFromEnum(sym);
  LifeState activeCells;
  activeCells.JoinWSymChain(activeRegion, symChain);
  std::vector<SymmetryTransform> waysToFlip;
  waysToFlip.push_back(Identity);
  switch(sym){
    // options are either of the two orthogonal reflections
    // each has the same effect.
    case StaticSymmetry::C2verticaleven:
      waysToFlip.push_back(ReflectAcrossY);
      break;
    case StaticSymmetry::C2horizontaleven:
      waysToFlip.push_back(ReflectAcrossX);
      break;
    case StaticSymmetry::C2:
      waysToFlip.push_back(ReflectAcrossX);
      waysToFlip.push_back(ReflectAcrossYeqX);
      waysToFlip.push_back(Rotate90);
      break;
    case StaticSymmetry::C2even:
      waysToFlip.push_back(ReflectAcrossXEven);
      waysToFlip.push_back(ReflectAcrossYeqX);
      waysToFlip.push_back(Rotate90Even);
      break;
    case StaticSymmetry::D2AcrossX:
      waysToFlip.push_back(ReflectAcrossY);
      waysToFlip.push_back(ReflectAcrossYEven);
      break;
    case StaticSymmetry::D2AcrossY:
      waysToFlip.push_back(ReflectAcrossX);
      waysToFlip.push_back(ReflectAcrossXEven);
      break;
    case StaticSymmetry::D2diagodd:
      waysToFlip.push_back(ReflectAcrossYeqNegX);
      waysToFlip.push_back(ReflectAcrossYeqNegXP1);
      break;
    case StaticSymmetry::D2negdiagodd:
    case StaticSymmetry::D4:
    case StaticSymmetry::D4even:
    case StaticSymmetry::C4:
    case StaticSymmetry::C4even:
      waysToFlip.push_back(ReflectAcrossYeqX);
      break;
  }
  
  int earliestGen = std::atoi(argv[5]);
  int latestGen = std::atoi(argv[6]);
  for ( auto trans : waysToFlip ){
    std::cout << "orfilter " << earliestGen << "-" << latestGen << " ";
    LifeState transformed = activeRegion;
    transformed.Transform(trans);
    std::cout << FilterForRegion(transformed) << std::endl;
  }
  // activeRegion.Print();

}
