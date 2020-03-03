#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <bits/stdc++.h> 
using namespace std;

int main() {
  string keyword;

  int vertexes[3][3];
  int backslashes;

  for(int i=0;i<3;i++) {
    cin>>keyword;
    backslashes = count(keyword.begin(), keyword.end(), '/');
    for(int j=0;j<3;j++) {
      if(keyword[0] != '/') {
        try {
        vertexes[i][j] = stoi(keyword)-1;
        } catch(...) {
          cout << "Invalid number " + keyword + "\n";
          exit(1);
        }
      }
      if(keyword.length() != 0) {
        keyword = keyword.substr(keyword.find("/") + 1);
      }
    }
    cout << "\n";
  }
  cout << "\n";
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      cout << vertexes[i][j] << " ";
    }
    cout << "\n";
  }
}