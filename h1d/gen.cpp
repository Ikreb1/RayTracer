#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

//to generate a 26x26 cube of sphere that all are next to eachother
//used for the imgToShow(was quite slow)
int main() {
  double fimm = 0.5;
  int len = 13;
  for(int i=-13;i<len;i++) {
    for(int j=-13;j<len;j++) {
      cout << "sphere " << fimm*i << " " << fimm*j << " 2.0 1.0\n";
    }
  }
}
