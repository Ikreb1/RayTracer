#include <iostream>
#include <string>
#include <fstream>
using namespace std;

	//makes a diamond thingy if height and width are multiple of 256
int main() {
	
	//takes in name of input and output file
	string inputFileName="";
	string outputFileName="";
	cin>>inputFileName;
	cin>>outputFileName;
	if(outputFileName=="") outputFileName="output";	
	
	//opens input file and sends out error if it doesn't work
	ifstream input;
	input.open("./" +inputFileName);
	if(!input) {
		cerr << "Unable to open file "+inputFileName+".txt\n";
		exit(1);
	}

	//initializing height and width
	//factor is used to make the image still pretty if you want some multiple of 256
	//if factor is 1 then the "correct" image only works if height and width are 256
	string item="";
	int factor = 4;
        int height;
        int width;
	int red = 0;
        int green = 0;
        int blue = 0;

	//Reads the input file and checks that it is in the correct format
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

	//Enlarges the image without ruining the formula	
	width = width*factor;
	height = height*factor;
	
	//creates an outpu file
	ofstream output;
	output.open(outputFileName+".ppm", ios::out | ios::trunc);
	
	if(output.is_open()) {
		//adds the required parameters
		output << "P3\n";
		output << width << "\n";
		output << height << "\n";
		output << 255 << "\n";

		//adds the colors and makes a pretty picture if the image is a factor of 256x256
		for(int i=0;i<width;i++) {
			for(int j=0;j<height;j++) {
				output << ( red + abs(i-width/2)/factor + abs(j-height/2)/factor)%256
				<< " " << (green)%256 
				<< " " << (blue + abs(i-width/2)/factor + abs(j-height/2)/factor)%256 
				<< " " <<"\n";
				red = red + 10;
				green = green + 2;
			}
		}
		//tries to close output file and throws error if not successful
		output.close();
	} else {
		cout << "Unable to open file"+ outputFileName +".ppm";
	}
	return 0;
}
