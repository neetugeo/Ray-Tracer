/*
CSCI 480
Assignment 3 Raytracer

Name: <Neetu George>
*/

#include <iostream>
#include <pic.h>
#include <windows.h>
#include <stdlib.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <stdio.h>
#include <string>
#include <math.h>
#include <list>
#include <vector>
#include <string>
#include <sstream>


#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;
bool firstPass = false;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode=MODE_JPEG;
int objectCounter = 0;
int counter = 0;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

std::list <struct Ray> rays;
std::vector <struct Ray> raysV;
float aspectRatio = WIDTH/HEIGHT;

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Direction
{
	float x;
	float y;
	float z;
};

struct Ray 
{
	struct Direction d;
	struct Direction normD;
	struct Direction origin;
	int coordX;
	int coordY;
};

struct Direction normalize(struct Direction d)
{

	float total = sqrt(pow(d.x,2) + pow(d.y,2) + pow(d.z,2));

	d.x = d.x/(total);
	d.y = d.y/(total);
	d.z = d.z/(total);

	return(d);
}

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);
float dot(Direction r1, Direction r2);
float intersectPlane (Ray r, Triangle t, Direction& n);


void saveScreenshot (char *filename)
{
		int i, j;
		Pic *in = NULL;

		if (filename == NULL)
		return;

		/* Allocate a picture buffer */
		in = pic_alloc(640, 480, 3, NULL);

		printf("File to save to: %s\n", filename);

		for (i=479; i>=0; i--) 
		{
			glReadPixels(0, 479-i, 640, 1, GL_RGB, GL_UNSIGNED_BYTE,
            &in->pix[i*in->nx*in->bpp]);
		}

		if (jpeg_write(filename, in))
			printf("File saved Successfully\n");
		else
			printf("Error in Saving\n");

		pic_free(in);	
}


void calculateRay()
{
  float trX, tlX, brX, blX;
  float trY, tlY, brY, blY;
  float xDistance, yDistance;
  float Z = -1.0;
  float aspectRatio = (float)WIDTH/(float)HEIGHT;

  //Top left corner
  tlX= -1.0*(aspectRatio)*tan((fov/2.0) * (3.14/180));
  tlY = tan((fov/2.0) * (3.14/180));
 
  //Y right corner
  trX = aspectRatio*(((fov/2.0) * (3.14/180)));
  trY = tan((fov/2.0) * (3.14/180));

  //Bottom Left corner
  blX = -1.0*(aspectRatio)*(tan((fov/2.0) * (3.14/180)));
  blY = -1.0*tan((fov/2.0) * (3.14/180));

  //Bottom Right corner
  brX = aspectRatio*tan((fov/2.0) * (3.14/180));
  brY = -1.0*tan((fov/2.0) * (3.14/180));

  xDistance = (trX*2.0f)/(float)WIDTH;
  yDistance = (trY*2.0f)/(float)HEIGHT;

  //Make a 2D array. Store each vector ray
  for (int i = 0; i < WIDTH; i++)
  {
	  for (int j = 0; j < HEIGHT; j++)
	  {
		  float offsetX = (float)i*xDistance;
		  float offsetY = (float)j*yDistance;

		  struct Ray tempRay;
		  struct Direction temp; 
		  struct Direction normTemp;
		  struct Direction origin;

		  origin.x = 0;
		  origin.y = 0;
		  origin.z = 0;

		  temp.x = blX + offsetX; //we start at the bottom left
		  temp.y = blY + offsetY;
		  temp.z = -1;

		  normTemp = normalize(temp);
		  tempRay.d = temp;
		  tempRay.normD = normTemp;
		  tempRay.origin = origin;
		  tempRay.coordX = i;
		  tempRay.coordY = j;

		  rays.push_back(tempRay);
		  raysV.push_back(tempRay);
		  
	  }
  }
}

Direction crossProduct(Direction a, Direction b)
{
	Direction result;
	result.x = (a.y * b.z) - (a.z * b.y);
	result.y = (a.z * b.x) - (a.x * b.z);
	result.z = (a.x * b.y) - (a.y * b.x);

	return (result);
}

void illuminateTriangle(Ray re, Triangle t, Light li, float alpha, float beta, float gamma)
{
	Ray n;
	n.d.x = (alpha * t.v[0].normal[0]) + (beta * t.v[1].normal[0]) + (gamma * t.v[2].normal[0]);
	n.d.y = (alpha * t.v[0].normal[1]) + (beta * t.v[1].normal[1]) + (gamma * t.v[2].normal[1]);
	n.d.z = (alpha * t.v[0].normal[2]) + (beta * t.v[1].normal[2]) + (gamma * t.v[2].normal[2]);

	float eSpec[3]; 
	eSpec[0]= (alpha * t.v[0].color_specular[0]) + (beta * t.v[1].color_specular[0]) + (gamma * t.v[2].color_specular[0]);
	eSpec[1]= (alpha * t.v[0].color_specular[1]) + (beta * t.v[1].color_specular[1]) + (gamma * t.v[2].color_specular[1]);
	eSpec[2]= (alpha * t.v[0].color_specular[2]) + (beta * t.v[1].color_specular[2]) + (gamma * t.v[2].color_specular[2]);

	float eDiff[3];
	eDiff[0]= (alpha * t.v[0].color_diffuse[0]) + (beta * t.v[1].color_diffuse[0]) + (gamma * t.v[2].color_diffuse[0]);
	eDiff[1]= (alpha * t.v[0].color_diffuse[1]) + (beta * t.v[1].color_diffuse[1]) + (gamma * t.v[2].color_diffuse[1]);
	eDiff[2]= (alpha * t.v[0].color_diffuse[2]) + (beta * t.v[1].color_diffuse[2]) + (gamma * t.v[2].color_diffuse[2]);

	float eShi;

	eShi= (alpha * t.v[0].shininess) + (beta * t.v[1].shininess) + (gamma * t.v[2].shininess);

	Ray l;
	l.d.x = li.position[0] - re.d.x;
	l.d.y = li.position[1] - re.d.y;
	l.d.z = li.position[2] - re.d.z;
	l.normD = normalize(l.d);



	Ray r;
	float c = 2.0f * (dot(l.normD,n.d));
	r.d.x = c*(n.d.x) - l.normD.x;
	r.d.y = c*(n.d.y) - l.normD.y;
	r.d.z = c*(n.d.z) - l.normD.z;

	struct Ray v;
	v.d.x = 0.0f - re.d.x;
	v.d.y = 0.0f - re.d.y;
	v.d.z = 0.0f - re.d.z;

	v.normD = normalize(v.d);
	float valueLD = dot(l.normD,n.d);
	float valueRV = dot(r.d,v.normD);
	
	if (valueLD < 0)
	{
		valueLD = 0;
	}
	if (valueRV < 0)
	{
		valueRV = 0;
	}

	float Ir = li.color[0] * (eDiff[0] * (valueLD) + eSpec[0] * pow((double)valueRV,(double)eShi));
	float Ig = li.color[1] * (eDiff[1] * (valueLD) + eSpec[1] * pow((double)valueRV,(double)eShi));
	float Ib = li.color[2] * (eDiff[2] * (valueLD) + eSpec[2] * pow((double)valueRV,(double)eShi));

	float Ar = ambient_light[0];
	float Ag = ambient_light[1];
	float Ab = ambient_light[2];

    Ir = Ir + Ar;
	Ig = Ig + Ag;
	Ib = Ib + Ab;

	if (Ir < 0)
	{
		Ir = 0;
	}
	if (Ig < 0)
	{
		Ig = 0;
	}
	if (Ib < 0)
	{
		Ib = 0;
	}

	if (Ir > 1)
	{
		Ir = 1;
	}
	if (Ig > 1)
	{
		Ig = 1;
	}
	if (Ib > 1)
	{
		Ib = 1;
	}

	plot_pixel(re.coordX,re.coordY,(Ir*255),(Ig*255),(Ib*255));
}

float triangleArea(Vertex a, Vertex b, Vertex c, Ray r, Direction n, Triangle t)
{
	if (n.x < 0)
		n.x = -1.0 * n.x;

	if (n.y < 0)
		n.y = -1.0 * n.y;

	if (n.z < 0)
		n.z = -1.0 * n.z;
 
	float largest = max(n.x,n.y);
	largest = max(n.z,largest);

	float area = 0.0f;
	
	if(largest == n.z)
	{
		float area = (0.5f) * ((b.position[0] - a.position[0]) * (c.position[1] - a.position[1]) - (c.position[0] - a.position[0]) * (b.position[1] - a.position[1]));
		return(area);
	}
	else if (largest == n.y)
	{
		float area = (0.5f) * ((b.position[0] - a.position[0]) * (c.position[2] - a.position[2]) - (c.position[0] - a.position[0]) * (b.position[2] - a.position[2]));
		return(area);
	}
	else if (largest == n.x)
	{
		float area = (0.5f) * ((b.position[1] - a.position[1]) * (c.position[2] - a.position[2]) - (c.position[1] - a.position[1]) * (b.position[2] - a.position[2]));
		return(area);
	}

	return(area);
	
}

void barycentric(Triangle t, float& t0, Ray r, float& alpha, float& beta, float& gamma)
{
	Direction n;
	t0 = intersectPlane(r, t, n);
	Vertex v;



	if (t0 > 0)
	{
		v.position[0] = (r.normD.x * t0) - r.origin.x;
		v.position[1] = (r.normD.y * t0) - r.origin.y;
		v.position[2] = (r.normD.z * t0) - r.origin.z;

		float TTA = triangleArea(t.v[0],t.v[1],t.v[2], r, n, t);
		float Cn12 = triangleArea(v, t.v[1], t.v[2], r, n, t);
		float C0n2 = triangleArea(t.v[0], v, t.v[2], r, n, t);
		float C01n = triangleArea(t.v[0], t.v[1], v, r, n, t);

		alpha = Cn12/TTA;
		beta = C0n2/TTA;
		gamma = C01n/TTA;
	}
}

float intersectPlane (struct Ray ray, Triangle t, Direction& normal)
{
	float t0 = -1.0f;

	Direction BA;
	Direction CA;

	BA.x = t.v[1].position[0] - t.v[0].position[0];
	BA.y = t.v[1].position[1] - t.v[0].position[1];
	BA.z = t.v[1].position[2] - t.v[0].position[2];

	CA.x = t.v[2].position[0] - t.v[0].position[0];
	CA.y = t.v[2].position[1] - t.v[0].position[1];
	CA.z = t.v[2].position[2] - t.v[0].position[2];

	normal = crossProduct(BA,CA);
	normal = normalize(normal);

	/*if (normal.z < 0)
	{
		normal.x = normal.x * (-1.0f);
		normal.y = normal.y * (-1.0f);
		normal.z = normal.z * (-1.0f);
	}*/
	
	float d = (-1.0f) * ((normal.x * t.v[0].position[0]) + (normal.y * t.v[0].position[1]) + (normal.z * t.v[0].position[2]));

	if (dot(normal,ray.normD) != 0)
	{
		t0 = ((-1.0f) * (dot(normal,ray.origin) + d)) / dot(normal,ray.normD); 
	}

	if (t0 <= 0)
	{
		return(-1.0f);
	}
	else if (t0 > 0) 
	{
		return(t0);
	}
	else
	{
		return(-1.0f);
	}

}


float intersectSphere(struct Ray ray, Sphere s)
{
	float a = 1;
	float b = 2.0f*(ray.normD.x * (ray.origin.x - s.position[0]) + ray.normD.y * (ray.origin.y - s.position[1]) + ray.normD.z * (ray.origin.z - s.position[2]));
	float c = pow((ray.origin.x - s.position[0]),2) + pow((ray.origin.y - s.position[1]),2) + pow((ray.origin.z - s.position[2]),2) - pow(s.radius,2);

	if ((pow(b,2) - 4.0f*c) >= 0)
	{
		float t0 = ((-1.0f)*b - sqrt((pow(b,2) - 4.0f*c)))/2.0f;
		float t1 = ((-1.0f)*b + sqrt((pow(b,2) - 4.0f*c)))/2.0f;

		if (t0 > 0 && t1 > 0)
		{
			return (min(t0,t1));
		}
		else
		{
			return(-1.0f);
		}
		
	}

	else
	{
		return(-1.0f);
	}
}

float dot(struct Direction r1, struct Direction r2)
{
	return ((r1.x * r2.x) + (r1.y * r2.y) + (r1.z * r2.z));
}

void illuminateSphere(struct Ray re, Sphere s, Light li)
{
	struct Ray n;
	n.d.x = re.d.x - s.position[0];
	n.d.y = re.d.y - s.position[1];
	n.d.z = re.d.z - s.position[2];
	n.d = normalize(n.d);


	struct Ray l;
	l.d.x = li.position[0] - re.d.x;
	l.d.y = li.position[1] - re.d.y;
	l.d.z = li.position[2] - re.d.z;
	l.normD = normalize(l.d);



	struct Ray r;
	float c = 2.0f * (dot(l.normD,n.d));
	r.d.x = c*(n.d.x) - l.normD.x;
	r.d.y = c*(n.d.y) - l.normD.y;
	r.d.z = c*(n.d.z) - l.normD.z;

	struct Ray v;
	v.d.x = 0.0f - re.d.x;
	v.d.y = 0.0f - re.d.y;
	v.d.z = 0.0f - re.d.z;

	v.normD = normalize(v.d);
	float valueLD = dot(l.normD,n.d);
	float valueRV = dot(r.d,v.normD);
	
	if (valueLD < 0)
	{
		valueLD = 0;
	}
	if (valueRV < 0)
	{
		valueRV = 0;
	}
	float Ir = li.color[0] * (s.color_diffuse[0] * (valueLD) + s.color_specular[0] * pow((double)valueRV,s.shininess));
	float Ig = li.color[1] * (s.color_diffuse[1] * (valueLD) + s.color_specular[1] * pow((double)valueRV,s.shininess));
	float Ib = li.color[2] * (s.color_diffuse[2] * (valueLD) + s.color_specular[2] * pow((double)valueRV,s.shininess));

	float Ar = ambient_light[0];
	float Ag = ambient_light[1];
	float Ab = ambient_light[2];

    Ir = Ir + Ar;
	Ig = Ig + Ag;
	Ib = Ib + Ab;

	if (Ir < 0)
	{
		Ir = 0;
	}
	if (Ig < 0)
	{
		Ig = 0;
	}
	if (Ib < 0)
	{
		Ib = 0;
	}

	if (Ir > 1)
	{
		Ir = 1;
	}
	if (Ig > 1)
	{
		Ig = 1;
	}
	if (Ib > 1)
	{
		Ib = 1;
	}

	plot_pixel(re.coordX,re.coordY,(Ir*255),(Ig*255),(Ib*255));
}

float distance(struct Ray r, float t0)
{
	float x = r.normD.x * t0;
	float y = r.normD.y * t0;
	float z = r.normD.z * t0;
	return (sqrt(pow(x,2) + pow(y,2) + pow(z,2)));
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;

  if (firstPass == false)
  {
	  calculateRay();
	  firstPass = true;
  }

  //simple output

  float tMin = 10000;
  
  glPointSize(2.0);  
  glBegin(GL_POINTS);

	for (int j = 0; j < raysV.size(); j++)
	 {
		for (int k = 0; k < num_lights; k++)
			  {
				float t0, t2;
				float alpha, beta, gamma;

				for(int i = 0; i < num_spheres; i++)
				{
					counter++;
					t0 = intersectSphere(raysV[j],spheres[i]);

					if (t0 != -1.0f)
					{
						struct Ray thisRay;
						thisRay.d.x = raysV[j].normD.x * t0;
						thisRay.d.y = raysV[j].normD.y * t0;
						thisRay.d.z = raysV[j].normD.z * t0;
						thisRay.normD = normalize(thisRay.d);
						thisRay.coordX = raysV[j].coordX;
						thisRay.coordY = raysV[j].coordY;

						struct Ray tempRay;
						tempRay.d.x = raysV[j].normD.x * t0;
						tempRay.d.y = raysV[j].normD.y * t0;
						tempRay.d.z = raysV[j].normD.z * t0;
						tempRay.normD = normalize(tempRay.d);


						tempRay.coordX = raysV[j].coordX;
						tempRay.coordY = raysV[j].coordY;
						tempRay.origin.x = thisRay.d.x;
						tempRay.origin.y = thisRay.d.y;
						tempRay.origin.z = thisRay.d.z;
						tempRay.d.x = lights[k].position[0] - thisRay.d.x;
						tempRay.d.y = lights[k].position[1] - thisRay.d.y;
						tempRay.d.z = lights[k].position[2] - thisRay.d.z;
						
						tempRay.normD = normalize(tempRay.d);

						for (int p = 0; p < num_spheres; p++)
						{
							float t1 = intersectSphere(tempRay, spheres[p]);
						
							if (t1 != -1)
							{
								if (p != i)
								{
									plot_pixel(raysV[j].coordX,raysV[j].coordY,0,0,0);
								}
							}
							else 
							{
								if (t0 < tMin)
								{
									tMin = t0;
									illuminateSphere(thisRay,spheres[i],lights[k]);
								}
							}
						}

						for (int p = 0; p < num_triangles; p++)
						{
							Direction n;
							float t1 = intersectPlane(tempRay, triangles[p], n);
						
							if (t1 != -1)
							{
								if (t0 < tMin)
								{
									tMin = t0;
									plot_pixel(raysV[j].coordX,raysV[j].coordY,0,0,0);
								}
								
							}
							else 
							{
								if (t0 < tMin)
								{
									tMin = t0;
									illuminateSphere(thisRay,spheres[i],lights[k]);
								}
							}
						} 
					 }
				  }
				for (int a = 0; a < num_triangles; a++)
				{
					counter++;
					barycentric(triangles[a],t2,raysV[j],alpha,beta,gamma);


					if (t2 != -1.0f)
					{
						struct Ray thisRay;
						thisRay.d.x = raysV[j].normD.x * t2;
						thisRay.d.y = raysV[j].normD.y * t2;
						thisRay.d.z = raysV[j].normD.z * t2;
						thisRay.normD = normalize(thisRay.d);
						thisRay.coordX = raysV[j].coordX;
						thisRay.coordY = raysV[j].coordY;

						struct Ray tempRay;
						tempRay.d.x = raysV[j].normD.x * t2;
						tempRay.d.y = raysV[j].normD.y * t2;
						tempRay.d.z = raysV[j].normD.z * t2;
						tempRay.normD = normalize(tempRay.d);


						tempRay.coordX = raysV[j].coordX;
						tempRay.coordY = raysV[j].coordY;
						tempRay.origin.x = thisRay.d.x;
						tempRay.origin.y = thisRay.d.y;
						tempRay.origin.z = thisRay.d.z;
						tempRay.d.x = lights[k].position[0] - thisRay.d.x;
						tempRay.d.y = lights[k].position[1] - thisRay.d.y;
						tempRay.d.z = lights[k].position[2] - thisRay.d.z;
						
						tempRay.normD = normalize(tempRay.d);

						float alp,bet,gam;
						for (int p = 0; p < num_spheres; p++)
						{
							float t1 = intersectSphere(tempRay,spheres[p]);
							if (t1 != -1)
							{
								if (t2 < tMin)
								{
									plot_pixel(raysV[j].coordX,raysV[j].coordY,0,0,0);
								}
							}
							else 
							{
								
								if (t2 < tMin)
								{
									
									
									if (alpha <= 0.0 && beta <= 0.0 && gamma <= 0.0 || alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0)
									{
										tMin = t2;
										illuminateTriangle(thisRay,triangles[a],lights[k],alpha,beta,gamma);
									}
								}
							}
						}
					
						for (int z = 0; z < num_triangles; z++)
						{
							//std::cout << "HERE" << std::endl;
							Direction n;
							float t1 = intersectPlane(tempRay,triangles[z], n);
							if (t1 != -1)
							{
								if (t2 < tMin)
								{
									plot_pixel(raysV[j].coordX,raysV[j].coordY,0,0,0);
								}
							}
							else 
							{
								
								if (t2 < tMin)
								{	
									if (alpha <= 0.0 && beta <= 0.0 && gamma <= 0.0 || alpha >= 0.0 && beta >= 0.0 && gamma >= 0.0)
									{
										tMin = t2;
										illuminateTriangle(thisRay,triangles[a],lights[k],alpha,beta,gamma);
									}
								}
							}
						}

					 }
				}
				tMin = 10000;
			  }
	 }
  glEnd();
  glFlush();
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;
  std::string* str = new std::string();
  *str = "001.jpg";

  char* filename = new char[str->length()+1];
  strcpy(filename,str->c_str());

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(stricmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(stricmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  triangles[num_triangles++] = t;
	}
      else if(stricmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
	  spheres[num_spheres++] = s;
	}
      else if(stricmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }

  objectCounter = num_spheres + num_triangles;

  return 0;
}

void display()
{

}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
	  mode = MODE_JPEG;
      draw_scene();
      if(mode == MODE_JPEG)
		save_jpg();
    }
  once=1;

}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
