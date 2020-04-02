#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <bits/stdc++.h> 
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

struct SimpleColorType {
  double r, g, b;
};

//stores a sphere cords and color
struct SphereType{
	double x, y, z;
	double r;
	//index of array where color is stored
	int m;
  //index of array where .ppm name is stored
  int t;
};

struct TextureCord {
  double u, v;
};

struct TriangleType{
  Point v1, v2, v3;
  Point vn1, vn2, vn3;
  TextureCord vt1, vt2, vt3;
  //n counts how many values are in triangle as in v, vn, vt
  int m, t, n;

  double alpha, beta, gamma;
};

struct TextureMapVector {
  int height, width;
  vector<vector<SimpleColorType> > v;
};
// vector<vetor<Color> > bydudde;

//vector is a resizable array to store color or sphere
vector<SphereType> sphereList;
//first color in list is background color
vector<ColorType> colorList;
//List of  lights
vector<Light> lightList;
//List of vertex points for triangles
vector<Point> vertexList;
//List of Triangles
vector<TriangleType> triangleList;
//List of vertex normal
vector<Point> vertexNormalList;
//list of .ppm file textures
vector<string> textureList;
//list of texture cordinates
vector<TextureCord> vertexTextureList;

vector<TextureMapVector> textureVector;


//function decleration
ColorType shadeRay(int id, int type,  Point intersect, Point dir, Point interpolatedN);
ColorType TraceRay(RayType ray);
double vectorLen(Point i);
Point unitVector(Point i);
Point vectorAddition(Point a, Point b);
Point vectorScalar(Point p, double x);
Point vectorFromPoints(Point p, Point q);
Point crossProduct(Point u, Point v);
Point unitVectorOfCrossProduct(RayType ray);
Point unitVectorOfCrossProduct(Point a, Point b);
Point makeAPoint(double x, double y, double z);
Point rayEquationVector(RayType ray, double t);
void printPoint(Point a, string s);
void printRay(RayType r, string s);
int isDouble(string s);
double dotProduct(Point a, Point b);
bool clampCheck(double x);
double isOutOfRange(double x);
bool isNearlyEqual(double x, double y);
int shadowRay(RayType shadow, int w, double d, int id, int type);
void printDouble(double x, string s);
ColorType toPolarCoorinateMapping(SphereType s, Point intersect);
ColorType toPolarCoorinateMapping(TriangleType t);


//returns color of sphere or background, type = 2 = triangle, type = 1 = sphere
ColorType shadeRay(int id, int type,  Point intersect, Point dir, Point interpolatedN) {
  //didn't hit any sphere so returns background color
  if(id == -1) {
    return colorList[0];
  }
  //extract the sphere we are calculation from the arrays
  
  ColorType color; 
  ColorType textureColor;
  SphereType sphere;
  TriangleType triangle;
  if(type == 1) {
    sphere = sphereList[id];
    color = colorList[sphere.m];
  } else {
    triangle = triangleList[id];
    color = colorList[triangle.m];
  }
  ColorType i;
  RayType shadow;
  Point l, v, n_v, h;
  Point p0, p1, p2, e1, e2;
  int shadowFlag;
  double ka = color.ka;
  double kd = color.kd;
  double ks = color.ks;
  int n = color.n;
  double d;
  double n_dot_l, n_dot_h;
  shadow.ox = intersect.x;
  shadow.oy = intersect.y;
  shadow.oz = intersect.z;
  //calc N vector outside loop for efficiency
  if(type == 1) {
    n_v = unitVector(makeAPoint((intersect.x - sphere.x)/sphere.r, (intersect.y - sphere.y)/sphere.r, (intersect.z - sphere.z)/sphere.r));
    if(sphere.t > -1) {
      
      textureColor = toPolarCoorinateMapping(sphere, intersect);
      color.odr = textureColor.odr;
      color.odg = textureColor.odg;
      color.odb = textureColor.odb;
    }
  } else {
    if(triangle.t > -1) {
      textureColor = toPolarCoorinateMapping(triangle);
      color.odr = textureColor.odr;
      color.odg = textureColor.odg;
      color.odb = textureColor.odb;
    }
    if(triangle.n != 1) {
      p0 = triangle.v1;
      p1 = triangle.v2;
      p2 = triangle.v3;
      e1 = vectorFromPoints(p0, p1);
      e2 = vectorFromPoints(p0, p2);
      n_v = crossProduct(e1, e2);
    } else {
      n_v = interpolatedN;
    }
  }
  v = unitVector(makeAPoint(-dir.x, -dir.y, -dir.z));
  //taking ka*od outside the loop since it only needs to be done once
  i.odr = ka*color.odr;
  i.odg = ka*color.odg;
  i.odb = ka*color.odb;
  for(Light light: lightList) {
    if(light.w == 0) {
      //change to L vector
      l = unitVector(makeAPoint(-light.x,-light.y,-light.z));
      //shadow direction for directional lights is the same
      shadow.dx = l.x;
      shadow.dy = l.y;
      shadow.dz = l.z;
      //test if our sphere is blocked by another sphere 
      shadowFlag = shadowRay(shadow, light.w, 0, id, type);
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
      shadowFlag = shadowRay(shadow, light.w, d, id, type);
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
    i.odr += shadowFlag*isOutOfRange(light.r*kd*color.odr + ks*color.osr);
    i.odg += shadowFlag*isOutOfRange(light.g*kd*color.odg + ks*color.osg);
    i.odb += shadowFlag*isOutOfRange(light.b*kd*color.odb + ks*color.osb);
    /*
    i.odr += light.r*isOutOfRange(isOutOfRange(kd*color.odr) + isOutOfRange(ks*color.osr));
    i.odg += light.g*isOutOfRange(isOutOfRange(kd*color.odg) + isOutOfRange(ks*color.osg));
    i.odb += light.b*isOutOfRange(isOutOfRange(kd*color.odb) + isOutOfRange(ks*color.osb));
    */
  }
  i.odr = isOutOfRange(i.odr);
  i.odg = isOutOfRange(i.odg);
  i.odb = isOutOfRange(i.odb);
	return i;
}

//finds sphere intersections from ray/sphere intersect points to light, so we can implement shadows
int shadowRay(RayType shadow, int w, double distance, int id, int type) {
  double A = 1;
	double B;
	double C;
	double discriminant;
	double t1, t2;
  int i = 0;
	//iterate through all spheres to find any sphere that is in the way of the light
	for(SphereType s: sphereList) {
    //ignores it self but should be changed to epsilon later
    if(!(i == id && type == 1)) {
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
          if(t1 <= distance && t1 >= 0)
            return 0;
          if(t2 <= distance && t2 >= 0)
            return 0;
        }
      }
    }
    i++;
	}

  Point e1, e2, p0, p1, p2, triN, p, ep, interpolatedN, vn1, vn2, vn3;
  double D, divisor, sanity, area; 
  double d11, d12, d22, dp1, dp2, d, dgamma, dbeta;
  double alpha, beta, gamma;
  //to check for intersections with triangles
  i = 0;
  for(TriangleType tri: triangleList) {
    if(!(type == 2 && i == id)) {
      t1 = -1;
      p0 = tri.v1;
      p1 = tri.v2;
      p2 = tri.v3;
      e1 = vectorFromPoints(p0, p1);
      e2 = vectorFromPoints(p0, p2);
      triN = crossProduct(e1, e2);
      //to see if we intersect the plane of the triangle
      D = -1.0 * (triN.x*p0.x + triN.y*p0.y + triN.z*p0.z);
      divisor = triN.x*shadow.dx + triN.y*shadow.dy + triN.z*shadow.dz;
      //if divisor == 0 we are done since there will be no intersection
      if(divisor != 0) {
        t1 = -1.0*((triN.x*shadow.ox + triN.y*shadow.oy + triN.z*shadow.oz + D)/divisor);
        //if t is less than zero we are also done since then the intersect  point is behind us
        if(!(t1 < 0)) {
          //does Xo + t*Xd to which is the ray/plane intersection point
          p = rayEquationVector(shadow, t1);
          
          area = 0.5 * vectorLen(triN);
          ep = vectorFromPoints(p0, p);
          d11 = dotProduct(e1, e1);
          d12 = dotProduct(e1, e2);
          d22 = dotProduct(e2, e2);
          dp1 = dotProduct(ep, e1);
          dp2 = dotProduct(ep, e2);
          d = d11*d22 - d12*d12;
          dbeta = d22*dp1 - d12*dp2;
          dgamma = d11*dp2 - d12*dp1;
          beta = dbeta / d;
          gamma = dgamma / d;
          alpha = 1.0 - (beta + gamma);
          
          if(!isNearlyEqual(alpha + beta + gamma, 1)) {
            //cout << "shadow: not right since alpha + beta + gamma = " << alpha + gamma + beta << "\n";
          }

          //checks if alpha, beta and gamma all are in the range 0-1 and does rounding errors aswell
          if(0 <= alpha && alpha <= 1 || isNearlyEqual(alpha, 0) || isNearlyEqual(alpha, 1)) {
            if(0 <= beta && beta <= 1 || isNearlyEqual(beta, 0) || isNearlyEqual(beta, 1)) {
              if(0 <= gamma && gamma <= 1 || isNearlyEqual(gamma, 0) || isNearlyEqual(gamma, 1)) {
                if(t1 >= 0) {
                  //if directional light we only need to know if distance is positive and since we already checked for negative
                  //we know there is a sphere in the way
                  if(w == 0) {
                    return 0;
                  }
                  //for point light we need to make sure that the sphere intersection isnt behind us and also isnt past the point light
                  if(t1 <= distance)
                    return 0;
                }
              }
            }
          }
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
  //make ray.dx be the correct thing that is the direction of the eye at x pixel point
	double len = sqrt( pow(ray.dx-ray.ox, 2) + pow(ray.dy-ray.oy, 2) + pow(ray.dz-ray.oz, 2));
	ray.dx = (ray.dx-ray.ox)/len;
	ray.dy = (ray.dy-ray.oy)/len;
	ray.dz = (ray.dz-ray.oz)/len;

	double A = 1;
	double B;
	double C;
	double discriminant;
	double t1, t2;
	double t = -1;
	double temp;
  int i = 0;
	int sphereId = -1;
  int triangleId = -1;
	//iterate through all spheres to find smallest t.
	for(SphereType s: sphereList) {
		temp = -1;
		//calc B and C for each sphere since A is always 1
		B = 2*(ray.dx*(ray.ox-s.x)+ray.dy*(ray.oy-s.y)+ray.dz*(ray.oz-s.z));
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
          sphereId = i;
        }
      }
		}
    i++;
	}

  Point e1, e2, p0, p1, p2, triN, p, ep, interpolatedN, vn1, vn2, vn3;
  double D, divisor, sanity, area; 
  double d11, d12, d22, dp1, dp2, d, dgamma, dbeta;
  double alpha, beta, gamma;
  //to check for intersections with triangles
  
  /*
  ray.ox = 0.0;
  ray.oy = 0.0;
  ray.oz = 0.0;
  ray.dx = 1.0/3;
  ray.dy = 2.0/3;
  ray.dz = 2.0/3;
  */
  i = 0;
  for(TriangleType tri: triangleList) {
    temp = -1;
    p0 = tri.v1;
    p1 = tri.v2;
    p2 = tri.v3;
    //printPoint(p0, "p0");
    //printPoint(p1, "p1");
    //printPoint(p2, "p2");
    e1 = vectorFromPoints(p0, p1);
    //printPoint(e1, "e1");
    e2 = vectorFromPoints(p0, p2);
    //printPoint(e2, "e2");
    triN = crossProduct(e1, e2);
    //printPoint(triN, "n");
    //to see if we intersect the plane of the triangle
    D = -1.0 * (triN.x*p0.x + triN.y*p0.y + triN.z*p0.z);
    //cout << "D = " << D << "\n";
    divisor = triN.x*ray.dx + triN.y*ray.dy + triN.z*ray.dz;
    //cout << "divisor = " << divisor << "\n";
    //if divisor == 0 we are done since there will be no intersection
    if(divisor != 0) {
      temp = -1.0*((triN.x*ray.ox + triN.y*ray.oy + triN.z*ray.oz + D)/divisor);
      //printDouble(temp, "t");
      //if t is less than zero we are also done since then the intersect  point is behind us
      if(!(temp < 0)) {
        //does Xo + t*Xd to which is the ray/plane intersection point
        p = rayEquationVector(ray, temp);
        //printPoint(p, "p");
        
        
        /*
        doesnt return 0 but should it?
        //sanity check incase I am a dumbass
        sanity = p.x + p.y + p.z + D;
        if(!isNearlyEqual(sanity, 0)) {
          cout << "something is wrong in the triangle intersection with plane equaton\n";
          cout << "might be a < epsilon problem since sanity = " << sanity << "\n";
          //exit(1);
        }
        */
        
        area = 0.5 * vectorLen(triN);
        ep = vectorFromPoints(p0, p);
        d11 = dotProduct(e1, e1);
        d12 = dotProduct(e1, e2);
        d22 = dotProduct(e2, e2);
        dp1 = dotProduct(ep, e1);
        dp2 = dotProduct(ep, e2);
        d = d11*d22 - d12*d12;
        dbeta = d22*dp1 - d12*dp2;
        dgamma = d11*dp2 - d12*dp1;
        beta = dbeta / d;
        gamma = dgamma / d;
        alpha = 1.0 - (beta + gamma);
        //printDouble(alpha, "alpha");
        //printDouble(beta, "beta");
        //printDouble(gamma, "gamma");
        if(!isNearlyEqual(alpha + beta + gamma, 1)) {
          cout << "not right since alpha + beta + gamma = " << alpha + gamma + beta << "\n";
        }

        //checks if alpha, beta and gamma all are in the range 0-1 and does rounding errors aswell
        if(0 <= alpha && alpha <= 1 || isNearlyEqual(alpha, 0) || isNearlyEqual(alpha, 1)) {
          if(0 <= beta && beta <= 1 || isNearlyEqual(beta, 0) || isNearlyEqual(beta, 1)) {
            if(0 <= gamma && gamma <= 1 || isNearlyEqual(gamma, 0) || isNearlyEqual(gamma, 1)) {
              if(t == -1 || t > temp) {
                t = temp;
                triangleId = i;
                if(tri.n == 1) {
                  //finding the interpolated N
                  vn1 = vectorScalar(tri.vn1, alpha);
                  vn2 = vectorScalar(tri.vn2, beta);
                  vn3 = vectorScalar(tri.vn3, gamma);
                  //adds the three vectors together and then returns unitvector of them
                  interpolatedN = unitVector(vectorAddition(vectorAddition(vn1, vn2), vn3));
                }
                if(tri.t > -1) {
                  
                  triangleList[i].beta = beta;
                  triangleList[i].alpha = alpha;
                  triangleList[i].gamma = gamma;
                }
              }
            } //else cout << "gamma is not inside 0-1 gamma = " << gamma << "\n";
          } //else cout << "beta is not inside 0-1 beta = " << beta << "\n";
        } //else cout << "alpha is not inside 0-1 alpha = " << alpha << "\n";
      }
    }
    i++;
  }
  
  intersectionPoint = makeAPoint(ray.ox + t*ray.dx, ray.oy + t*ray.dy, ray.oz + t*ray.dz);
  if(triangleId != -1) {
    return shadeRay(triangleId, 2, intersectionPoint, dir, interpolatedN);
  } else {
    //intersection point = point(pxview.o_lambda) + t*vector(dir)
    dir = makeAPoint(ray.dx, ray.dy, ray.dz);
    //now we find the color of the sphere and if any lights are illumination it and return that color
    return shadeRay(sphereId, 1, intersectionPoint, dir, dir);
  }
}

//changes the spheres XYZ center to polar coordinates
ColorType toPolarCoorinateMapping(SphereType s, Point intersect) {
  //inside acos and atan needs to be in range -1, +1
  //returns -Pi to Pi rad
  double phi =  acos((intersect.z - s.z) / s.r);
  //returns -PI/2, PI/2 rad

  double theta = atan2((intersect.y - s.y), (intersect.x - s.x));
  //acos returns values in -Pi to Pi rad not 0 - PI see slides 92
  
  //v is now 0-1
  double v = phi / M_PI;
  //is u correct? should be range 0-1
  double u = (theta + M_PI)/(2*M_PI);
  int height = (textureVector[s.t].height) - 1;
  int width = (textureVector[s.t].width) - 1;
  int placementv = (int) (v*height);
  int placementu = (int) (u*width);
  //cout << placementu << " " << placementv << "\n";
  
  SimpleColorType sc = textureVector[s.t].v.at(placementv).at(placementu);
  //cout << u << ": u " << v << ":v \n";
  //cout << sc.r << " " << sc.b << " " << sc.g << "\n";
  ColorType c;
  c.odr = sc.r;
  c.odg = sc.g;
  c.odb = sc.b;
  //cout << c.odr << " " << c.odb << " " << c.odg << "\n";
  return c;
}

ColorType toPolarCoorinateMapping(TriangleType t) {
  TextureCord temp;
  //cout << t.alpha << " " << t.vt1.u << "\n";
  double u = t.alpha * t.vt1.u + t.beta * t.vt2.u + t.gamma * t.vt3.u;
  double v = t.alpha * t.vt1.v + t.beta * t.vt2.v + t.gamma * t.vt3.v;
  
  //cout << u << " " << v << "\n";

  int height = textureVector[t.t].height - 1;
  int width = textureVector[t.t].width - 1;
  int placementv = (int) (v*height);
  int placementu = (int) (u*width);
  //cout << placementu << " " << placementv << "\n";
  SimpleColorType sc = textureVector[t.t].v.at(placementu).at(placementv);
  
  ColorType c;
  c.odr = sc.r;
  c.odg = sc.g;
  c.odb = sc.b;
  
  return c;
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

void printDouble(double x, string s) {
  cout << s << " = " << x << "\n";
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

//does vector addtion of 2 vectors
Point vectorAddition(Point a, Point b) {
  Point temp;
  temp.x = a.x + b.x;
  temp.y = a.y + b.y;
  temp.z = a.z + b.z;
  return temp;
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

//checks if double number is pretty much equal since sometimes they are off by 10^-20 or something
bool isNearlyEqual(double x, double y) {
  double epsilon = 0.000001;
  //cout << "abs(" << x << " - " << y << ") <= " << epsilon << "\n";
  //cout << "abs(x - y) = " << abs(x - y) << "\n";
  return abs(x - y) <= epsilon;
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

//multiplies a vector with a constant
Point vectorScalar(Point p, double x) {
  p.x = p.x * x;
  p.y = p.y * x;
  p.z = p.z * x;
  return p;
}

//gets a vector from 2 points
Point vectorFromPoints(Point a, Point b) {
  return makeAPoint(b.x-a.x, b.y-a.y, b.z-a.z);
}

//cross product UxV = (uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx)
Point crossProduct(Point u, Point v) {
	Point temp;
	temp.x = u.y*v.z-u.z*v.y;
	temp.y = u.z*v.x-u.x*v.z;
	temp.z = u.x*v.y-u.y*v.x;	return temp;
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

//calculates ray = ray.ox + t*ray.dx for xyz
Point rayEquationVector(RayType ray, double t) {
  Point p;
  p.x = ray.ox + t*ray.dx;
  p.y = ray.oy + t*ray.dy;
  p.z = ray.oz + t*ray.dz;
  return p;
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
  TriangleType triangle;
  Point vertex;
  Point vertexNormal;
  TextureCord vertexTexture;
  Light light;
	int matColorId = 0;
  int textureId = -1;
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
			if(matColorId == 0) {
        cout << "You need to define a material color before sphere" << "\n";
        exit(1);
      }
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

      sphere.t = textureId;
			sphere.m = matColorId;
			sphereList.push_back(sphere);
		} else if(keyword=="v"){
      input>>keyword;
			try {
			vertex.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			vertex.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			vertex.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
      vertexList.push_back(vertex); 
    } else if(keyword=="vn"){
      input>>keyword;
			try {
			vertexNormal.x = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			vertexNormal.y = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			vertexNormal.z = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}
      vertexNormalList.push_back(vertexNormal);
      } else if(keyword=="vt"){
      input>>keyword;
			try {
			vertexTexture.u = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

			input>>keyword;
			try {
			vertexTexture.v = stod(keyword);
			} catch(...) {
				cout << "Invalid number " + keyword + "\n";
				exit(1);
			}

      vertexTextureList.push_back(vertexTexture); 
    } else if(keyword=="f"){
      //vertexes[i][j] i=input number, j=v/vt/vn
      int vertexes[3][3];
      int backslashes;
      int normal = 0;

      input>>keyword;
      backslashes = count(keyword.begin(), keyword.end(), '/');
      if(backslashes == 2) {
        for(int i=0;i<3;i++) {
          if(i != 0) input>>keyword;
          for(int j=0;j<3;j++) {
            if(keyword[0] != '/') {
              try {
              vertexes[i][j] = stoi(keyword)-1;
              if(stoi(keyword) == 0) {
              cout << "vertex id start from 1 not 0\n";
              exit(1);
              }
              } catch(...) {
                cout << "Invalid number " + keyword + "\n";
                exit(1);
              }
            } else {
              vertexes[i][j] = -1;
            }
            if(keyword.length() != 0) {
              keyword = keyword.substr(keyword.find("/") + 1);
            }
          }
        }
      } else if(backslashes == 1) {
        for(int i=0;i<3;i++) {
          if(i != 0) input>>keyword;
          for(int j=0;j<3;j++) {
            if(j == 2) vertexes[i][j] = -1;
            else if(keyword[0] != '/') {
              try {
              vertexes[i][j] = stoi(keyword)-1;
              if(stoi(keyword) == 0) {
              cout << "vertex id start from 1 not 0\n";
              exit(1);
              }
              } catch(...) {
                cout << "Invalid number " + keyword + "\n";
                exit(1);
              }
            } else {
              vertexes[i][j] = -1;
            }
            if(keyword.length() != 0) {
              keyword = keyword.substr(keyword.find("/") + 1);
            }
          }
        }
      } else if(backslashes == 0) {
        for(int i=0;i<3;i++) {
          if(i != 0) input>>keyword;
          if(keyword[0] != '/') {
            try {
            vertexes[i][0] = stoi(keyword)-1;
            if(stoi(keyword) == 0) {
              cout << "vertex id start from 1 not 0\n";
              exit(1);
            }
            } catch(...) {
              cout << "Invalid number " + keyword + "\n";
              exit(1);
            }
          } else {
            vertexes[i][0] = -1;
          }
          if(keyword.length() != 0) {
            keyword = keyword.substr(keyword.find("/") + 1);
          }
          vertexes[i][1] = -1;
          vertexes[i][2] = -1;
        }
      } else {
        cout << "too many '/' in argument: " << keyword << "\n";
        exit(1);
      }

      /*
      for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
          cout << vertexes[i][j] << " ";
        }
        cout << "\n";
      }

      /*
      for(int i=0;i<3;i++) {
        if(i!=1 && (vertexes[i][0] == vertexes[i][1] || vertexes[i][1] == vertexes[i][2] || vertexes[i][2] == vertexes[i][0])) {
          cout << i << "\n";
          cout << "A triangle can't have 2 of the same vertex\n";
          exit(1);
        }
      }
      */
      
      //checks if value is out of range
      int size = vertexList.size()-1;
      for(int j=0;j<3;j++) {
        if(vertexes[j][0] < -1 || vertexes[j][0] > size) {
          cout << "v" << 0 << " and value " << j << " is out of range\n";
          exit(1);
        }
      }

      size = vertexTextureList.size()-1;
      for(int j=0;j<3;j++) {
        if(vertexes[j][1] < -1 || vertexes[j][1] > size) {
          cout << "v" << 1 << " and value " << j << " is out of range\n";
          cout << vertexes[j][1] << "\n";
          exit(1);
        }
      }

      size = vertexNormalList.size()-1;
      for(int j=0;j<3;j++) {
        if(vertexes[j][2] < -1 || vertexes[j][2] > size) {
          cout << "v" << 2 << " and value " << j << " is out of range\n";
          exit(1);
        }
      }
			
      triangle.v1 = vertexList[vertexes[0][0]];
			triangle.v2 = vertexList[vertexes[1][0]];
			triangle.v3 = vertexList[vertexes[2][0]];
      if(vertexes[0][2] != -1) {
        triangle.vn1 = vertexNormalList[vertexes[0][2]];
        triangle.vn2 = vertexNormalList[vertexes[1][2]];
        triangle.vn3 = vertexNormalList[vertexes[2][2]];
        normal++;
      }
      if(vertexes[0][1] != -1) {
        triangle.vt1 = vertexTextureList[vertexes[0][1]];
        triangle.vt2 = vertexTextureList[vertexes[1][1]];
        triangle.vt3 = vertexTextureList[vertexes[2][1]];
      }
      triangle.n = normal;
      triangle.t = textureId;
      triangle.m = matColorId;
      triangleList.push_back(triangle);
    } else if(keyword == "texture") {
      input>>keyword;
      textureList.push_back(keyword);
      textureId++;
    } else {
      cerr << "unknown argument \"" + keyword + "\" in inputfile\n";
      exit(1);
		}    
	}
  input.close();
	//checks for parallel of up and view_dir
	if(up.x==viewdir.x&&up.y==viewdir.y&&up.z==viewdir.z) {
		cout << "viewdir and up can not be parallel\n";
		exit(1);
	}


  for(int i=0;i<textureList.size();i++) {
    //read in texture input files, temp only for one
    ifstream texture;
    TextureMapVector textureMapVector;
    texture.open("./" + textureList[i]);
    if(!texture) {
      cerr << "Unable to open file "+textureList[0]+"\n";
      exit(1);
    }

    texture>>keyword;
    if(keyword != "P3") {
      cout << "texture file isn't formated as a P3 file";
      exit(1);
    }

    texture>>keyword;
    textureMapVector.width = stoi(keyword);

    texture>>keyword;
    textureMapVector.height = stoi(keyword);
    texture>>keyword;
    if(stoi(keyword) != 255) {
      cout << textureList[0] << ".ppm file should be a 255 value format it is: " << keyword << "\n";
      exit(1);
    }

    //SimpleColorType textureColor[textureHeight][textureWidth];
    //vector<SimpleColorType **> bruh;
    //bruh.push_back(textureColor **);
    
    
    SimpleColorType simpleColor;
    vector<SimpleColorType> colorVector;
    for(int i=0;i<textureMapVector.height;i++) {
      for(int j=0;j<textureMapVector.width;j++) {
          texture>>keyword;
          simpleColor.r = stod(keyword)/255;
          texture>>keyword;
          simpleColor.g = stod(keyword)/255;
          texture>>keyword;
          simpleColor.b = stod(keyword)/255;
          colorVector.push_back(simpleColor);
      }
      textureMapVector.v.push_back(colorVector);
      colorVector.clear();
    }
    textureVector.push_back(textureMapVector);
    textureMapVector.v.clear();
    texture.close();
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
