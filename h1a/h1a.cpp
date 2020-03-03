#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
using namespace std;

//ease of storage of vector/point
struct Point{
	double x,y,z;
};

//the ray from eye to a point in view if i understand correctly
struct RayType{
	//origin
	double ox, oy, oz;
	//destination
	double dx, dy, dz;
};

//stores the rgb
struct ColorType{
	double red, green, blue;
};

//stores a sphere cords and color
struct SphereType{
	double x, y, z;
	double r;
	//pointer to array where color is stored
	int m;
};

//vector is a resizable array to store color or sphere
vector<SphereType> sphereList;
//first color in list is background color
vector<ColorType> colorList;


//function decleration
ColorType shadeRay(int colorId);
ColorType TraceRay(RayType ray, int i, int j, double *array, int width);
double vectorLen(Point i);
Point unitVector(Point i);
Point crossProduct(Point u, Point v);
Point unitVectorOfCrossProduct(RayType ray);
Point unitVectorOfCrossProduct(Point a, Point b);
Point makeAPoint(double x, double y, double z);
void printPoint(Point a, string s);
int isDouble(string s);


//returns color of sphere or background
ColorType shadeRay(int colorId) {
	/* for assignment 1a, parameters should include:
	the coordinates of the point where the shading
	will be computed, and an identifier that allows
	access to the material properties of the object

	compute the color at the point; presently this
	consists of just returning the sphere’s material
	color. Later on the illumination calculations
	will become more complicated, and this is where
	the recursion will take place when your program
	is extended from ray casting to ray tracing */
	//Return computed_color
		return colorList[colorId];
}

//calcs if there is an intersection. doesn't work i think
ColorType TraceRay(RayType ray, int i, int j, double *array, int width) {
	//How to define each ray
	//view_orign = eye_xyz
	//P = ray.dx
	//dir = Pxview_origin //unit vector
	//then we find ray
	//if t=0 we are at the same place, - is behind and t>(Pxview_orign.len) means
	//we are beyond view window
	
	Point dir;
	double len = sqrt( pow(ray.dx-ray.ox, 2) + pow(ray.dy-ray.oy, 2) + pow(ray.dz-ray.oz, 2));
	dir.x = (ray.dx-ray.ox)/len;
	dir.y = (ray.dy-ray.oy)/len;
	dir.z = (ray.dz-ray.oz)/len;

	/*
	cout << "x: " << ray.dx << "\n";
	cout << "y: " << ray.dy << "\n";
	cout << "z: " << ray.dz << "\n";
	*/
	//printPoint(dir, "dir");
	//ray = view_origin + t*(dir)
	//Point rayray = makeAPoint(ray.ox+t*dir.x, ray.oy+t*dir.y, ray.oz+t*dir.z);

	//intersection of sphere and ray
	//t = unknown
	//ray equation = (ray.ox, ray.oy, ray.oz)+t*(ray.dx, ray.dy, ray.dz)
	//x=ray.x+t*ray.dx y=... etc.
	//sphere equation is:
	//(x-x_c)²+(y-y_c)²+(z-z_c)² = r²
	//put ray eq in sphere eq
	//then we want to solve for t
	//1. A*t²+B*t+C=0
	//A=(ray.dx²+ray.dy²+ray.dz²)=1//since they are unit vector
	//B=2*(ray.dx*(ray.x-sphere.x)+ray.dy*(ray.y-sphere.y)+ray.dz*(ray.z-sphere.z))
	//C=(ray.x-sphere.x)²+(ray.z-sphere.z)²+(ray.z-sphere.z)² - sphere.r
	//solve eq 1 with (-B(+/-)sqrt(B²-4*A*C))/2*A
	//whats inside the sqrt is called discriminant
	//if inside sqrt is positive then they intersect and we can see the sphere
	//if inside sqrt = 0 then we are grazing the sphere
	//if inside sqrt is neg we missed the sphere
	double A = 1;
	double B;
	double C;
	double discriminant;
	double t1, t2;
	array[i*width+j]=-1;
	double temp;
	int colorId = 0;

	//iterate through all spheres to find smallest t. doesnt work
	for(SphereType s: sphereList) {
		temp = -1;
		//calc B and C for each sphere since A is always 1
		B = 2*(dir.x*(ray.ox-s.x)+dir.y*(ray.oy-s.y)+dir.z*(ray.oz-s.z));
		C = pow((ray.ox-s.x), 2.0) + pow((ray.oz-s.z), 2.0) + pow((ray.oz-s.z), 2.0) - s.r;
		//calc discriminant of quadratic formula
		discriminant = B*B-4*A*C;
		//cout << discriminant << "\n";
		//if the ray and sphere intersect we try to find the color of it, if the sphere is the closest one.
		if(discriminant >= 0) {
			t1 = (-B + sqrt(discriminant))/2*A;
			t2 = (-B - sqrt(discriminant))/2*A;
			//cout << "t1: " << t1 << "\n";
			//cout << "t2: " << t2 << "\n";
			if(t1 >= 0) {
				if(t2 < 0) {temp = t1;} 
				else if(t1 < t2) {temp = t2;} 
				else {temp = t1;}
				//cout << "found one " << t1 << " " << t2 << "\n";
			} else if(t2 >= 0) {
				if(t1 < 0) {temp = t2;}
				else if(t2 < t1) {temp = t1;}
				else {temp = t2;}
				//cout << "found one " << t1 << " " << t2 << "\n";
			}
			if(temp != -1 && array[i*width+j] < temp) {
				array[i*width+j] = temp;
				colorId = s.m;
			}
		} //else cout << "no hit\n";
	}
	//cout << "\n";
	return shadeRay(colorId);
}
/*
“ray” specifies the incident ray (origin and direction); 
	for assignment 1a it is the ray from the eye through a point 
	in the view window. please be forewarned that this parameter list will
	evolve; later on we will need to pass in more information 
	
ColorType Trace_Ray( RayType ray ) {
	for each object in the scene, check for a ray/object intersection; 
	keep track of the closest intersection point and of the identity
	of the object hit at that point
	call Shade_Ray() to determine the color 
	return( return_color );
}
*/

//shorthand to debug
void printPoint(Point a, string s) {
	cout << s << "\n";
	cout << "x: " << a.x << "\n";
	cout << "y: " << a.y << "\n";
	cout << "z: " << a.z << "\n";
}

//vector length calc
double vectorLen(Point i) {
	return sqrt(i.x*i.x+i.y*i.y+i.z*i.z);
}

//makes a vector into a unit vector
Point unitVector(Point i) {
	Point temp;
	double len = vectorLen(i);
	temp.x = i.x/len;
	temp.y = i.y/len;
	temp.z = i.z/len;
	return temp;
}

//cross product UxV = (uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx)
Point crossProduct(Point u, Point v) {
	Point temp;
	temp.x = u.y*v.z-u.z*v.y;
	temp.y = u.z*v.x-u.x*v.z;
	temp.z = u.x*v.y-u.y*v.x;
	return temp;
}

//cross product of d and o in ray and then converts it to unit vector
Point unitVectorOfCrossProduct(RayType ray) {
	Point a = makeAPoint(ray.dx, ray.dy, ray.dz);
	Point b = makeAPoint(ray.ox, ray.oy, ray.oz);
	return unitVectorOfCrossProduct(a, b);
}

//cross product of axb and then converts outcome to unit vector
Point unitVectorOfCrossProduct(Point a, Point b) {
	return unitVector(crossProduct(a, b));
}

//short hand to fill in a point
Point makeAPoint(double x, double y, double z) {
	Point temp;
	temp.x = x;
	temp.y = y;
	temp.z = z;
	return temp;
}



// complier command: gcc h1a.cpp -o h1a -lm -lstdc++
int main(int argc, char *argv[]) {
	//to do
	/* read scene description from input file */
	/* initialize pixel array for output image */
	/* perform necessary preliminary calculations */ 
	
	/* for each pixel in the image array: */
		/* call Trace_Ray() with appropriate parameters */
		/* use the value returned by Trace_Ray() to update the pixel color */
	
	/* write the final image to an output file */
    
	//takes in name of input and output file

	//varoiable decleration
	string inputFileName="";
	string outputFileName="";
    string keyword;
	RayType currentPixelPlace;
	ColorType backgroundColor;
	ColorType currentPixelColor;
	ColorType materialColor;
	int matColorId = 0;
	SphereType sphere;
	Point viewdir;
	//view_origin
	Point eye;
	//up dir cords
	Point up;
	//w = -view_dir
	Point n;
	//field of view
	double fov;
	//image size
	int height;
	int width;
	double aspectRatio;

	Point u;
	Point v;
	//size of view
	double viewWidth;
	double viewHeight;
	Point ul, ur, ll, lr;
	Point viewCenter;
	Point deltach;
	Point deltacv;
	Point deltah;
	Point deltav;
	double t;
	double d;

	
	//just notes i took when reading slides ignore it
	//vetor len = sqrt(x*x+y*y+z*z)
	//cross product UxV = (uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx)
	//unit vector = (x/len, y/len, z/len)
	//remember the right hand rule for cross product
	//u = ( view_dir x up_dir )/u_vector_len
	//v = (u x view_dir)/v_vector len
	//dont forget about the right order, view then up

	//corners of viewing window
	//aspect_ratio = width/height
	//d can be anything?
	//w = 2*d*tan(fov/2)
	//h = aspect ratio of image
	//h = w / (width/height)
	//ul = view_origin + d*(view_dir/len) - (w/2)*u + h/2*v
	//ur = view_origin + d*(view_dir/len) + (w/2)*u + h/2*v
	//ll = view_origin + d*(view_dir/len) - (w/2)*u - h/2*v
	//lr = view_origin + d*(view_dir/len) + (w/2)*u - h/2*v

	//center of view window
	//p = eye + d*(unit_vector_view_dir)

	//pixel mapping
	//delta_C_h = (ur -ul)/(2*imsize_width)
	//delta_v_h = (ll -ul)/(2*imsize_height)
	//delta_h = (ur-ul)/imsize_width
	//delta_v = (ll-ul)/imsize_height
	//to determine color at each pixel(i,j)
	//ul+(i)*delta_h+(j)*delta_v+delta_c_h+delta_c_v

	//How to define each ray
	//view_orign is eye_xyz?
	//dir = Pxview_origin //unit vector
	//then we find ray
	//ray = view_origin + t*(dir)
	//if t=0 we are at the same place, - is behind and t>(Pxview_orign.len) means
	//we are beyond view window

	//intersection of sphere and ray
	//t = unknown
	//ray equation = (ray.x, ray.y, ray.z)+t*(ray.dx, ray.dy, ray.dz)
	//x=ray.x+t*ray.dx y=... etc.
	//sphere equation is:
	//(x-x_c)²+(y-y_c)²+(z-z_c)² = r²
	//put ray eq in sphere eq
	//then we want to solve for t
	//1. A*t²+B*t+C=0
	//A=(ray.dx²+ray.dy²+ray.dz²)=1//since they are unit vector
	//B=2*(ray.dx*(ray.x-sphere.x)+ray.dy*(ray.y-sphere.y)+ray.dz*(ray.z-sphere.z))
	//C=(ray.x-sphere.x)²+(ray.z-sphere.z)²+(ray.z-sphere.z)² - sphere.r
	//solve eq 1 with (-B(+/-)sqrt(B²-4*A*C))/2*A
	//whats inside the sqrt is called discriminant
	//if inside sqrt is positive then they intersect and we can see the sphere
	//if inside sqrt =0 then we are grazing the sphere
	//if inside sqrt is neg we missed the sphere

	//later we care about the t and we want the solution that is positive and smallest

	cin>>inputFileName;
	//cin>>outputFileName;  //outputfile name not used anymore?
	//if(outputFileName=="")
	outputFileName="output";	
	
	//opens input file and sends out error if it doesn't work
	ifstream input;
	input.open("./" +inputFileName);
	if(!input) {
		cerr << "Unable to open file "+inputFileName+"\n";
		exit(1);
	}

    //reads inputs and initializes them correctly and checks for errors
    while(input>>keyword) {
		if(keyword=="eye"){
			input>>keyword;
			try {
				eye.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			eye.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			eye.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
		} else if(keyword=="imsize") {
			input>>width;
			if(width<1){
				cerr << "width needs to be a non negative number\n";
				exit(1);
			}	
			input>>height;
			if(height<1){
				cerr << "height needs to be a non negative number\n";
			}
		} else if(keyword=="viewdir"){
			input>>keyword;
			try {
			viewdir.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			viewdir.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			viewdir.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

		} else if(keyword=="updir"){
			input>>keyword;
			try {
			up.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
			
			input>>keyword;
			try {
			up.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
			
			input>>keyword;
			try {
			up.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

		} else if(keyword=="hfov"){
			input>>keyword;
			try {
			fov = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			if(fov > 180) {
				cout << "hfov can not be greater then 180.0\n";
				exit(1);
			}
			//degrees to radiant
			fov=fov*M_PI/180;

		} else if(keyword=="bkgcolor"){
			input>>keyword;
			
			try {
			backgroundColor.red = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			backgroundColor.green = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			backgroundColor.blue = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			colorList.push_back(backgroundColor);
		} else if(keyword=="mtlcolor"){
			input>>keyword;
			try {
			materialColor.red = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			materialColor.green = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			materialColor.blue = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			colorList.push_back(materialColor);
			matColorId++;
		} else if(keyword=="sphere"){
			input>>keyword;
			try {
			sphere.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			sphere.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			sphere.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			sphere.r = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			sphere.m = matColorId;
			sphereList.push_back(sphere);
		} else {
			cerr << "unknown argument \"" + keyword + "\" in inputfile\n";
			exit(1);
		}
	}

	//checks for parallel of up and view_dir
	if(up.x==viewdir.x&&up.y==viewdir.y&&up.z==viewdir.z) {
		cout << "viewdir and up can not be parallel\n";
		exit(1);
	}
	ofstream output;
	output.open("output.ppm", ios::out | ios::trunc);
	if(!output.is_open()) {
		cout << "Unable to open file output.ppm";
		exit(1);
	}
	//adds the required parameters
	output << "P3\n";
	output << width << "\n";
	output << height << "\n";
	output << 255 << "\n";

	
	//remember the right hand rule for cross product
	//u = ( view_dir x up_dir )/u_vector_len
	//v = (u x view_dir)/v_vector len
	//calc w, u and v
	n = makeAPoint(-viewdir.x, -viewdir.y, -viewdir.z);
	u = unitVectorOfCrossProduct(viewdir, up);
	v = unitVectorOfCrossProduct(u, viewdir);
	
	//corners of viewing window
	//aspect_ratio = width/height
	//d can be anything?
	//w = 2*d*tan(fov/2)
	//h = aspect ratio of image
	//h = w / (width/height)
	//ul = view_origin + d*(view_dir/len) - (w/2)*u + h/2*v
	//ur = view_origin + d*(view_dir/len) + (w/2)*u + h/2*v
	//ll = view_origin + d*(view_dir/len) - (w/2)*u - h/2*v
	//lr = view_origin + d*(view_dir/len) + (w/2)*u - h/2*v
	
	double array[height*width];
	//calc all the constants that we need
	aspectRatio = (double)width/height;
	d = 2;///what?
	viewWidth = 2*d*tan(fov/2);
	//cout << "fov: " << fov << "\n";
	//cout << "tan(fov/2): " << tan(fov/2) << "\n";
	viewHeight = viewWidth/aspectRatio;
	
	double viewdirlen = vectorLen(viewdir);
					//forward dir              //left right //up down
	ul.x = eye.x + d*(viewdir.x/viewdirlen) - viewWidth/2*u.x + viewHeight/2*v.x;
	ul.y = eye.y + d*(viewdir.y/viewdirlen) - viewWidth/2*u.y + viewHeight/2*v.y;
	ul.z = eye.z + d*(viewdir.z/viewdirlen) - viewWidth/2*u.z + viewHeight/2*v.z;
	
	ur.x = eye.x + d*(viewdir.x/viewdirlen) + viewWidth/2*u.x + viewHeight/2*v.x;
	ur.y = eye.y + d*(viewdir.y/viewdirlen) + viewWidth/2*u.y + viewHeight/2*v.y;
	ur.z = eye.z + d*(viewdir.z/viewdirlen) + viewWidth/2*u.z + viewHeight/2*v.z;
	
	ll.x = eye.x + d*(viewdir.x/viewdirlen) - viewWidth/2*u.x - viewHeight/2*v.x;
	ll.y = eye.y + d*(viewdir.y/viewdirlen) - viewWidth/2*u.y - viewHeight/2*v.y;
	ll.z = eye.z + d*(viewdir.z/viewdirlen) - viewWidth/2*u.z - viewHeight/2*v.z;
	
	lr.x = eye.x + d*(viewdir.x/viewdirlen) + viewWidth/2*u.x - viewHeight/2*v.x;
	lr.y = eye.y + d*(viewdir.y/viewdirlen) + viewWidth/2*u.y - viewHeight/2*v.y;
	lr.z = eye.z + d*(viewdir.z/viewdirlen) + viewWidth/2*u.z - viewHeight/2*v.z;

	//center of view window
	//p = eye + d*(unit_vector_view_dir)
	viewCenter.x = eye.x + d*(viewdir.x/viewdirlen);
	viewCenter.y = eye.y + d*(viewdir.y/viewdirlen);
	viewCenter.z = eye.z + d*(viewdir.z/viewdirlen);

	//pixel mapping
	//delta_C_h = (ur-ul)/(2*imsize_width)
	//delta_v_h = (ll-ul)/(2*imsize_height)
	//delta_h = (ur-ul)/imsize_width
	//delta_v = (ll-ul)/imsize_height
	//calc the pixel mapping
	deltach.x = (ur.x - ul.x)/(2*width);
	deltach.y = (ur.y - ul.y)/(2*width);
	deltach.z = (ur.z - ul.z)/(2*width);

	deltacv.x = (ll.x - ul.x)/(2*height);
	deltacv.y = (ll.y - ul.y)/(2*height);
	deltacv.z = (ll.z - ul.z)/(2*height);
	
	deltah.x = (ur.x - ul.x)/width;
	deltah.y = (ur.y - ul.y)/width;
	deltah.z = (ur.z - ul.z)/width;
	
	deltav.x = (ll.x - ul.x)/height;
	deltav.y = (ll.y - ul.y)/height;
	deltav.z = (ll.z - ul.z)/height;
	
	/*
	printPoint(u, "u");
	printPoint(v, "v");
	cout << "viewwidth " << viewWidth << "\n";
	cout << "viewHeight " << viewHeight << "\n";
	cout << "viewdirlen " << viewdirlen << "\n";
	printPoint(ul, "ul");
	printPoint(ur, "ur");
	printPoint(ll, "ll");
	printPoint(lr, "lr");
	printPoint(deltach, "delta_c_h");
	printPoint(deltacv, "delta_c_v");
	printPoint(deltah, "delta_h");
	printPoint(deltav, "delta_v");
	cout << "\n";
	*/

	//to determine color at each pixel(i,j)
	//ul+(i)*delta_h+(j)*delta_v+delta_c_h+delta_c_v

	//How to define each ray
	//view_orign is eye_xyz 
	//dir = Pxview_origin //unit vector
	//then we find ray
	//ray = view_origin + t*(dir)
	//if t=0 we are at the same place, - is behind and t>(Pxview_orign.len) means
	//we are beyond view window
	 
	//since each ray has eye as origin we can do that before loop
	currentPixelPlace.ox = eye.x;
	currentPixelPlace.oy = eye.y;
	currentPixelPlace.oz = eye.z;
	//for each pixel we calculate points place in view and then run traceray to find what color is needed and
	//then fill in the output image with the correct rgb values
	for(int h=0;h<height;h++) {
		for(int j=0;j<width;j++) {
			currentPixelPlace.dx = ul.x + h*deltah.x + j*deltav.x + deltach.x+deltacv.x;
			currentPixelPlace.dy = ul.y + h*deltah.y + j*deltav.y + deltach.y+deltacv.y;
			currentPixelPlace.dz = ul.z + h*deltah.z + j*deltav.z + deltach.z+deltacv.z;
			currentPixelColor = TraceRay(currentPixelPlace, h, j, array, width);
			output << (int)(currentPixelColor.red*255) << " "
			<< (int)(currentPixelColor.green*255) << " "
			<< (int)(currentPixelColor.blue*255) << "\n";
		}
	}

	output.close();

    return 0;
}