#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
using namespace std;

int main() {
  double BOX_SIZE = 18.21533455;
  int num_windows = 20;

  double step_length = BOX_SIZE / (2*num_windows + 1);
  for (int i=0; i<num_windows; i++)
	cout << setprecision(10) << 2*i*step_length << "\t\t" << (2*i + 3)*step_length << endl;
  return 0;
}
