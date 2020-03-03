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

//light struct
struct Light {
  double x,y,z;
  int w;
  double r, g, b;
};

//A ray's origin and its direction
struct RayType{
	//origin
	double ox, oy, oz;
	//destination
	double dx, dy, dz;
};

//stores the rgb
struct ColorType{
	double odr, odg, odb;
  double osr, osg, osb;
  double ka, kd, ks;
  int n;
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
//List of  lights
vector<Light> lightList;

//function decleration
ColorType shadeRay(int sId, Point intersect, Point dir);
ColorType TraceRay(RayType ray);
double vectorLen(Point i);
Point unitVector(Point i);
Point crossProduct(Point u, Point v);
Point unitVectorOfCrossProduct(RayType ray);
Point unitVectorOfCrossProduct(Point a, Point b);
Point makeAPoint(double x, double y, double z);
void printPoint(Point a, string s);
void printRay(RayType r, string s);
int isDouble(string s);
double dotProduct(Point a, Point b);
bool clampCheck(double x);
double isOutOfRange(double x);
int shadowRay(RayType shadow, int w, double d, int sId);

//returns color of sphere or background
ColorType shadeRay(int sId, Point intersect, Point dir) {
  //didn't hit any sphere so returns background color
  if(sId == -1) {
    return colorList[0];
  }
  //extract the sphere we are calculation from the arrays
  SphereType sphere = sphereList[sId];
  ColorType sphereColor = colorList[sphere.m];
  ColorType i;
  RayType shadow;
  Point l, v, n_v, h;
  int shadowFlag;
  double ka = sphereColor.ka;
  double kd = sphereColor.kd;
  double ks = sphereColor.ks;
  int n = sphereColor.n;
  double d;
  double n_dot_l, n_dot_h;
  shadow.ox = intersect.x;
  shadow.oy = intersect.y;
  shadow.oz = intersect.z;
  //calc N vector outside loop for efficiency
  
  n_v = unitVector(makeAPoint((intersect.x - sphere.x)/sphere.r, (intersect.y - sphere.y)/sphere.r, (intersect.z - sphere.z)/sphere.r));
  v = unitVector(makeAPoint(-dir.x, -dir.y, -dir.z));
  //taking ka*od outside the loop since it only needs to be done once
  i.odr = ka*sphereColor.odr;
  i.odg = ka*sphereColor.odg;
  i.odb = ka*sphereColor.odb;
  for(Light light: lightList) {
    if(light.w == 0) {
      //change to L vector
      l = unitVector(makeAPoint(-light.x,-light.y,-light.z));
      //shadow direction for directional lights is the same
      shadow.dx = l.x;
      shadow.dy = l.y;
      shadow.dz = l.z;
      //test if our sphere is blocked by another sphere 
      shadowFlag = shadowRay(shadow, light.w, 0, sId);
    } else {
      //for point lights take unitvector of light-intersection points
      l = makeAPoint(light.x-intersect.x, light.y-intersect.y, light.z-intersect.z);
      d = vectorLen(l);
      l = unitVector(l);
      //shadow.dlambda = L
      shadow.dx = l.x;
      shadow.dy = l.y;
      shadow.dz = l.z;
      //test if our sphere is blocked by another sphere
      shadowFlag = shadowRay(shadow, light.w, d, sId);
    }
    //H for dot products
    h = unitVector(makeAPoint(l.x + v.x, l.y + v.y, l.z + v.z));
    
    //calc dot products
    n_dot_h = dotProduct(n_v, h);
    n_dot_l = dotProduct(n_v, l);
    //checks if the number is less than 0
    if(n_dot_l < 0) n_dot_l = 0;
    if(n_dot_h < 0) n_dot_h = 0;
    //multiplying the dot product and the k constants so we don't need to clac them for each color
    kd = kd*n_dot_l;
    ks = ks*pow(n_dot_h, n);
    //should be all the clamping we need
    i.odr += shadowFlag*isOutOfRange(light.r*kd*sphereColor.odr + ks*sphereColor.osr);
    i.odg += shadowFlag*isOutOfRange(light.g*kd*sphereColor.odg + ks*sphereColor.osg);
    i.odb += shadowFlag*isOutOfRange(light.b*kd*sphereColor.odb + ks*sphereColor.osb);
    /*
    i.odr += light.r*isOutOfRange(isOutOfRange(kd*sphereColor.odr) + isOutOfRange(ks*sphereColor.osr));
    i.odg += light.g*isOutOfRange(isOutOfRange(kd*sphereColor.odg) + isOutOfRange(ks*sphereColor.osg));
    i.odb += light.b*isOutOfRange(isOutOfRange(kd*sphereColor.odb) + isOutOfRange(ks*sphereColor.osb));
    */
  }
  i.odr = isOutOfRange(i.odr);
  i.odg = isOutOfRange(i.odg);
  i.odb = isOutOfRange(i.odb);
	return i;
}

//finds sphere intersections from ray/sphere intersect points to light, so we can implement shadows
int shadowRay(RayType shadow, int w, double d, int sId) {
  double A = 1;
	double B;
	double C;
	double discriminant;
	double t1, t2;
  int i = 0;
	//iterate through all spheres to find any sphere that is in the way of the light
	for(SphereType s: sphereList) {
    //ignores it self but should be changed to epsilon later
    if(i != sId) {
      //calc B and C for each sphere since A is always 1
      B = 2*(shadow.dx*(shadow.ox-s.x)+shadow.dy*(shadow.oy-s.y)+shadow.dz*(shadow.oz-s.z));
      C = pow((shadow.ox-s.x), 2.0) + pow((shadow.oy-s.y), 2.0) + pow((shadow.oz-s.z), 2.0) - pow(s.r, 2.0);
      //calc discriminant of quadratic formula
      discriminant = B*B-4*A*C;
      //if the ray and sphere intersect we know the origin is in shadow if t1 or t2 is positive
      if(discriminant >= 0) {
        t1 = (-B + sqrt(discriminant))/2*A;
        t2 = (-B - sqrt(discriminant))/2*A;
        if(t1 >= 0 && t2 >= 0) {
          //if directional light we only need to know if distance is positive and since we already checked for negative
          //we know there is a sphere in the way
          if(w == 0) {
            return 0;
          }
          //for point light we need to make sure that the sphere intersection isnt behind us and also isnt past the point light
          if(t1 <= d && t1 >= 0)
            return 0;
          if(t2 <= d && t2 >= 0)
            return 0;
        }
      }
    }
    i++;
	}
  //if no sphere in the way we can ignore shadow flag and return 1
  return 1;
}


//calcs if there is an intersection. doesn't work i think
ColorType TraceRay(RayType ray) {
	Point intersectionPoint;
	Point dir;
   
	double len = sqrt( pow(ray.dx-ray.ox, 2) + pow(ray.dy-ray.oy, 2) + pow(ray.dz-ray.oz, 2));
	dir.x = (ray.dx-ray.ox)/len;
	dir.y = (ray.dy-ray.oy)/len;
	dir.z = (ray.dz-ray.oz)/len;

	double A = 1;
	double B;
	double C;
	double discriminant;
	double t1, t2;
	double t = -1;
	double temp;
  int i = 0;
	int id = -1;
	//iterate through all spheres to find smallest t.
	for(SphereType s: sphereList) {
		temp = -1;
		//calc B and C for each sphere since A is always 1
		B = 2*(dir.x*(ray.ox-s.x)+dir.y*(ray.oy-s.y)+dir.z*(ray.oz-s.z));
		C = pow((ray.ox-s.x), 2.0) + pow((ray.oy-s.y), 2.0) + pow((ray.oz-s.z), 2.0) - pow(s.r, 2.0);
		//calc discriminant of quadratic formula
		discriminant = B*B-4*A*C;
    
		//if the ray and sphere intersect we try to find the color of it, if the sphere is the closest one.
		if(discriminant >= 0) {
      t1 = (-B + sqrt(discriminant))/2*A;
			t2 = (-B - sqrt(discriminant))/2*A;
      //to find if there is any positive t1 or t2 that is smaller than our current t
      if(t1 >= 0) {
				if(t2 < 0) {temp = t1;} 
				else if(t1 < t2) {temp = t1;} 
				else {temp = t2;}
			} else if(t2 >= 0) {
				if(t1 < 0) {temp = t2;}
				else if(t2 < t1) {temp = t2;}
				else {temp = t1;}
			}
      if(temp != -1) {
        if (t == -1 || t > temp) {
          t = temp;
          id = i;
        }
      }
		}
    i++;
	}
  //intersection point = point(pxview.o_lambda) + t*vector(dir)
  intersectionPoint = makeAPoint(ray.ox + t*dir.x, ray.oy + t*dir.y, ray.oz + t*dir.z);
  
  //now we find the color of the sphere and if any lights are illumination it and return that color
	return shadeRay(id, intersectionPoint, dir);
}

//shorthand to debug
void printPoint(Point a, string s) {
	cout << s << "\n";
	cout << fixed << "x: " << a.x << "\n";
	cout << fixed << "y: " << a.y << "\n";
	cout << fixed << "z: " << a.z << "\n";
  cout << "\n";
}

void printRay(RayType r, string s) {
  cout << s << "\n";
	cout << fixed << "ox: " << r.ox << "\n";
	cout << fixed << "oy: " << r.oy << "\n";
	cout << fixed << "oz: " << r.oz << "\n";
	cout << fixed << "dx: " << r.dx << "\n";
	cout << fixed << "dy: " << r.dy << "\n";
	cout << fixed << "dz: " << r.dz << "\n";
  cout << "\n";
}

//vector length calc
double vectorLen(Point i) {
	return sqrt(i.x*i.x+i.y*i.y+i.z*i.z);
}

//computes dot Product of two vectors
double dotProduct(Point a, Point b) {
  //cout << a.x << "*" << b.x << "+" << a.y << "*" << b.y << "+" << a.z << "*" << b.z << "\n";
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

//check that number is between 0 and 1
bool clampCheck(double x) {
  return x < 0 || x > 1;
}

//checks if number is out of range of 0-1 and clamps it to 0 or 1 if so
double isOutOfRange(double x) {
  if(x < 0) return 0.0;
  else if(x > 1) return 1.0;
  else return x;
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
	//takes in name of input and output file

	//variable decleration
	string inputFileName="";
	string outputFileName="";
  string keyword;
	RayType currentPixelPlace;
	ColorType backgroundColor;
	ColorType currentPixelColor;
	ColorType materialColor;
  Light light;
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
	Point ul, ur, ll;
	Point viewCenter;
	Point deltach;
	Point deltacv;
	Point deltah;
	Point deltav;
	double t;
	double d;

  Point tempPoint;

	cin>>inputFileName;
	outputFileName = inputFileName;	
  //outputFileName = "output"
	
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
			backgroundColor.odr = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			backgroundColor.odg = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			backgroundColor.odb = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      //checking if colors are in range of 0-1
      if(clampCheck(backgroundColor.odr) || clampCheck(backgroundColor.odg) || clampCheck(backgroundColor.odb)) {
        cout << "material color values (rgb) need to be less than or equal 1 and greater than or equal 0" << "\n";
        exit(1);
      }

			colorList.push_back(backgroundColor);
		} else if(keyword=="mtlcolor"){
			input>>keyword;
			try {
			materialColor.odr = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			materialColor.odg = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			materialColor.odb = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.osr = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.osg = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.osb = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.ka = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.kd = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.ks = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			materialColor.n = stoi(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
      if(materialColor.n < 1) {
        cout << "n can not be less than 1" << "\n";
        exit(1);
      }

      //checking if colors are in range of 0-1
      if(clampCheck(materialColor.odr) || clampCheck(materialColor.odg) || clampCheck(materialColor.odb)) {
        cout << "material color values (rgb) need to be less than or equal 1 and greater than or equal 0" << "\n";
        exit(1);
      }

      if(clampCheck(materialColor.osr) || clampCheck(materialColor.osg) || clampCheck(materialColor.osb)) {
        cout << "material color values (rgb) need to be less than or equal 1 and greater than or equal 0" << "\n";
        exit(1);
      }

			colorList.push_back(materialColor);
			matColorId++;

      //cout << "matC odr" << materialColor.odr << " odg:" << materialColor.odg << " odb: " << materialColor.odb << " osr:" << materialColor.osr << " osg:" << materialColor.osg 
      //<< " osb:" << materialColor.osb << " ka:" << materialColor.ka << " kd: " << materialColor.kd << " ks: " << materialColor.ks << " n: " << materialColor.n << "\n";
    } else if(keyword=="light"){
			input>>keyword;
			try {
			light.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			light.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			light.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			light.w = stoi(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      if(light.w != 1 && light.w != 0) {
        cout << "light w value needs to be 1 or 0 " << light.w << "\n";
        exit(1);
      }

      input>>keyword;
			try {
			light.r = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
      

      input>>keyword;
			try {
			light.g = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      input>>keyword;
			try {
			light.b = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      //checking if colors are in range of 0-1
      if(clampCheck(light.r) || clampCheck(light.g) || clampCheck(light.b)) {
        cout << "light color values (rgb) need to be less than or equal 1 and greater than or equal 0" << "\n";
        exit(1);
      }

      //cout << "light base x: " << light.x << " y: " << light.y << " z: " << light.z << " w: " << light.w << " r: " << light.r << " g: " << light.g << " b: " << light.b << "\n";

      if(light.w == 1) {
			  lightList.push_back(light);
        //cout << "hello\n";
      } else {
        lightList.push_back(light);
      }
			//cout << "light after x: " << light.x << " y: " << light.y << " z: " << light.z << " w: " << light.w << " r: " << light.r << " g: " << light.g << " b: " << light.b << "\n"; 
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
  //this is the v but V is -v
	v = unitVectorOfCrossProduct(u, viewdir);
	
	//calc all the constants that we need
	aspectRatio = (double)width/height;
	d = 2;
	viewWidth = 2*d*tan(fov/2);
	viewHeight = viewWidth/aspectRatio;
	double viewdirlen = vectorLen(viewdir);

  //center of view window
	//p = eye + d*(unit_vector_view_dir)
	viewCenter.x = eye.x + d*(viewdir.x/viewdirlen);
	viewCenter.y = eye.y + d*(viewdir.y/viewdirlen);
	viewCenter.z = eye.z + d*(viewdir.z/viewdirlen);

					//forward dir              //left right //up down
	ul.x = viewCenter.x - viewWidth/2*u.x + viewHeight/2*v.x;
	ul.y = viewCenter.y - viewWidth/2*u.y + viewHeight/2*v.y;
	ul.z = viewCenter.z - viewWidth/2*u.z + viewHeight/2*v.z;
	
	ur.x = viewCenter.x + viewWidth/2*u.x + viewHeight/2*v.x;
	ur.y = viewCenter.y + viewWidth/2*u.y + viewHeight/2*v.y;
	ur.z = viewCenter.z + viewWidth/2*u.z + viewHeight/2*v.z;
	
	ll.x = viewCenter.x - viewWidth/2*u.x - viewHeight/2*v.x;
	ll.y = viewCenter.y - viewWidth/2*u.y - viewHeight/2*v.y;
	ll.z = viewCenter.z - viewWidth/2*u.z - viewHeight/2*v.z;

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
	
	//to determine color at each pixel(i,j)
	//ul+(i)*delta_h+(j)*delta_v+delta_c_h+delta_c_v

	//since each ray has eye as origin we can do that before loop
	currentPixelPlace.ox = eye.x;
	currentPixelPlace.oy = eye.y;
	currentPixelPlace.oz = eye.z;
	//for each pixel we calculate points place in view and then run traceray to find what color is needed and
	//then fill in the output image with the correct rgb values
 	for(int h=0;h<height;h++) {
		for(int j=0;j<width;j++) {
			currentPixelPlace.dx = ul.x + j*deltah.x + h*deltav.x + deltach.x+deltacv.x;
			currentPixelPlace.dy = ul.y + j*deltah.y + h*deltav.y + deltach.y+deltacv.y;
			currentPixelPlace.dz = ul.z + j*deltah.z + h*deltav.z + deltach.z+deltacv.z;
			currentPixelColor = TraceRay(currentPixelPlace);
      output << (int)(currentPixelColor.odr*255) << " "
			<< (int)(currentPixelColor.odg*255) << " "
			<< (int)(currentPixelColor.odb*255) << "\n";
		}
	}

	output.close();

  return 0;
}