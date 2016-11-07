/// Title:    Testing different STL features, usage
/// Author:   Kyle M. Lang
/// Created:  2016-NOV-06
/// Modified: 2016-NOV-06

#include <iostream>
#include <algorithm>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

vector<int> getObsRows(int, vector<int>*);

int main() {

  vector<int> tmp1 = {1, 2, 3, 4, 5};
  vector<int> tmp2 = {1, 2, 3};
  vector<int> missArray[] = {tmp1, tmp2};
  vector<int> out = getObsRows(20, missArray);
  vector<int>::iterator it;

  cout << "Here are your observed rows, buddy:" << endl;
  for(it = out.begin(); it != out.end(); ++it)
    cout << *it << " ";
  cout << endl;

  return 0;
}


vector<int> getObsRows(int nObs, vector<int>* missArray) {
  vector<int> allRows(nObs);
  std::iota(allRows.begin(), allRows.end(), 1);
  vector<int> missRows = missArray[0];
  
  vector<int>::iterator it;
  vector<int> obsRows(max(missRows.size(), allRows.size()));
  
  it = set_difference(allRows.begin(),
		      allRows.end(),
		      missRows.begin(),
		      missRows.end(),
		      obsRows.begin());
  
  obsRows.resize(it - obsRows.begin());
  
  return out;
}

 
