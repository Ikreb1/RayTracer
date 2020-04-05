#include <iostream>
using namespace std;

int main() {
  cout << "#include <iostream>\n";
  cout << "using namespace std;\n\n";
  cout << "int main() {\n";
  
  for(int i=0;i<1001;i++) {
    cout << "cout << \"this is line number \" << " << i << " << \"\\n\";\n";
  }
  cout << "{\n";
}