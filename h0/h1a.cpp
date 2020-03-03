#include <iostream>
#include <string>
#include <fstream>
using namespace std;

typedef struct {
	float x, y, z;
	float dx, dy, dz;
} RayType;

typedef struct {

} ColorType;

typedef struct {
	float x, y, y;
	float r;
	int m;
} Spheretype;

ColorType Shade_Ray() {

}


int main() {
    	
	//takes in name of input and output file
	string inputFileName="";
	string outputFileName="";
    string keyword;
    double 

	cin>>inputFileName;
	cin>>outputFileName;
	if(outputFileName=="") outputFileName="output";	
	
	//opens input file and sends out error if it doesn't work
	ifstream input;
	input.open("./" +inputFileName);
	if(!input) {
		cerr << "Unable to open file "+inputFileName+"\n";
		exit(1);
	}

    //reads inputs and initializes them correctly and checks for errors
    while(input>>item) {
		if(item=="imsize") {

			
			input>>width;
			if(width<1){
				cerr << "width needs to be a non negative number\n";
				exit(1);
			}	
			input>>height;
			if(height<1){
				cerr << "height needs to be a non negative number\n";
			}
		} else { 
			cerr << "imsize is not the first word in inputfile\n"; 
			exit(1);
		}
		if(input>>item){
			cerr << "Too many arguments in inputfile\n";
			exit(1);
		}
	}

    return 0;
}