#include <TROOT.h>
#include <TMath.h>
#include <iostream>
#include <vector>

#include "inputs.h"


template<class T>
struct CompareIndicesByVectorValues {
  std::vector<T> *_values;
public:
  CompareIndicesByVectorValues(std::vector<T> *v) : _values(v)
  {}

  bool operator()(const int &i, const int &j)
  { return _values->at(i) > _values->at(j); }
};

template<class T>
void prnArr(const char *msg, const std::vector<T> &v)
{
  std::cout << msg << " [" << v.size() << "]: ";
  for (unsigned int i=0; i<v.size(); i++) {
    std::cout << " " << v[i];
  }
  std::cout << "\n";
}

template<class T>
void prnArr(const char *msg, const std::vector<T> &v,
		  const std::vector<int> &idx)
{
  std::cout << msg << " [" << v.size() << "]: ";
  for (unsigned int i=0; i<v.size(); i++) {
    std::cout << " " << v[idx[i]];
  }
  std::cout << "\n";
}

void testSort()
{
  std::vector<double> arr;
  arr.push_back(2.1);
  arr.push_back(3.);
  arr.push_back(1.);
  arr.push_back(4);

  std::vector<int> idx;
  for (unsigned int i=0; i<arr.size(); i++) idx.push_back(i);

  CompareIndicesByVectorValues<double> arrCmp(&arr);
  std::sort(idx.begin(),idx.end(),arrCmp);
  prnArr("arr",arr);
  prnArr("idx",idx);
  prnArr("arr sorted",arr,idx);


  std::cout << "test embeded code\n";
  std::vector<int> idx2=getSortDescendingIdx(arr);
  prnArr("arr",arr);
  prnArr("idx2",idx2);
  prnArr("arr sorted 2",arr,idx2);
}
