#include <GL/glut.h>
#include "bevgrafmath2017.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include "functions.cpp"
#include "optim/optim.hpp"

#include "IL/il.h"
#include <IL/ilu.h>
#include "LTexture.h"

const int winWidth = 800, winHeight = 600;
std::vector<vec2> points, optpoints;

int dragged = -1;
int pointRadius = 4;
int uDist = winWidth/4.0;
bool optimized = false;

//File loaded texture
LTexture gLoadedTexture;

void init() {
	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, winWidth, 0.0, winHeight);

	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable( GL_TEXTURE_2D );

	ilInit();
    ilClearColour( 255, 255, 255, 000 );
}


void dragPoint(std::vector<vec2> points, int sensitivity, int mouseX, int mouseY) {
	for (int i = 4; i < points.size(); i++) {
		if (dist(points[i], vec2(mouseX, winHeight - mouseY)) < sensitivity)
			dragged = i;
	}
}

double getYofLine(vec2 p1, vec2 p2, double x) {
	vec2 norm = vec2(p1.y - p2.y, -(p1.x -p2.x));

	return (norm.x * p1.x + norm.y * p1.y - norm.x * x) / norm.y;

}

void update(int val) {
	glutPostRedisplay();
	glutTimerFunc(2, update, 0);
}

void createU(int uDist) {
	vec2 o = points[0],
		x = points[1],
		y = points[2],
		z = points[3];
	
	//uX = o + uDist * normalize(x - o);
	//uY = o + uDist * normalize(y - o);
	//uZ = o + uDist * normalize(z - o);
	points.push_back(o + uDist * normalize(x - o));
	points.push_back(o + uDist * normalize(y - o));
	points.push_back(o + uDist * normalize(z - o));
	
}

void processMouse(int button, int action, int xMouse, int yMouse) {
	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN) {
		dragPoint(points, 10, xMouse, yMouse);
		if (dragged == -1 && points.size() < 4) {
			points.push_back(vec2(xMouse, winHeight - yMouse));
			if (points.size() == 4)
				createU(uDist);
		}
			
	}
		
	if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
		dragged = -1;
}



bool trueProj(vec2 o, vec2 x, vec2 y, vec2 z, vec2 uX, vec2 uY, vec2 uZ, double sen) {



	double alpha = acos(dot(uY - uX, uZ - uX)/ (length(uY - uX) * length(uZ - uX)));
	double beta = acos(dot(uX - uY, uZ - uY) / (length(uX - uY) * length(uZ - uY)));
	double gamma = acos(dot(uY - uZ, uX - uZ) / (length(uY - uZ) * length(uX - uZ)));

	double e = length(points[0] - points[1]);
	double g = length(points[0] - points[2]);
	double i = length(points[0] - points[3]);

	double f = length(points[1] - points[4]);
	double h = length(points[2] - points[5]);
	double j = length(points[3] - points[6]);

	if (std::abs(pow(e / f, 2) / pow(g / h, 2) - tan(alpha) / tan(beta)) > sen)
		return false;

	if (std::abs(pow(g / h, 2) / pow(i / j, 2) - tan(beta) / tan(gamma)) > sen)
		return false;

	if (std::abs(pow(e / f, 2) / pow(i / j, 2) - tan(alpha) / tan(gamma)) > sen)
		return false;


	return true;
}

void processMouseActiveMotion(int xMouse, int yMouse) {
	if (dragged >= 0) {
		points[dragged].x = xMouse;
		points[dragged].y = getYofLine(points[0], points[dragged-3], xMouse);

		std::cout << trueProj(points[0], points[1], points[2], points[3], points[4], points[5], points[6], 0.1) << std::endl;
	}
}

void displayPoints() {

	glColor3f(1.0, 0.0, 0.0);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < points.size(); i++) {
		glVertex2f(points[i].x, points[i].y);
	}
	glEnd();

}

void displayOptPoints(){
    glColor3f(0.0, 1.0, 0.0);
	glPointSize(12.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < optpoints.size(); i++) {
		glVertex2f(optpoints[i].x, optpoints[i].y);
	}
	glEnd();
        
}

void calcOptPoints(double t1, double t2, double t3){
	optpoints.clear();

    vec2 p1 = t1 * points[0] + (1-t1) * points[4];
    vec2 p2 = t2 * points[0] + (1-t2) * points[5];
    vec2 p3 = t3 * points[0] + (1-t3) * points[6];
    
    optpoints.push_back(p1);
    optpoints.push_back(p2);
    optpoints.push_back(p3);
    
    optimized = true;
}

void displayLines() {
	glColor3d(0.0, 0.0, 0.0);
	glLineWidth(2.0);
	
	for (int i = 4; i < 7; i++) {
		glBegin(GL_LINES);
		glVertex2d(points[0].x, points[0].y);
		glVertex2d(points[i].x, points[i].y);
		glEnd();
	}

	glBegin(GL_LINE_LOOP);
	glVertex2d(points[5].x, points[5].y);
	glVertex2d(points[1].x, points[1].y);
	glVertex2d(points[6].x, points[6].y);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex2d(points[5].x, points[5].y);
	glVertex2d(points[3].x, points[3].y);
	glVertex2d(points[4].x, points[4].y);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex2d(points[6].x, points[6].y);
	glVertex2d(points[2].x, points[2].y);
	glVertex2d(points[4].x, points[4].y);
	glEnd();
	
}

void displayOptLines() {
	glColor3d(0.0, 0.0, 1.0);
	glLineWidth(3.0);
	
	for (int i = 4; i < 7; i++) {
		glBegin(GL_LINES);
		glVertex2d(points[0].x, points[0].y);
		glVertex2d(points[i].x, points[i].y);
		glEnd();
	}

	for (int i = 4; i < 7; i++) 
	{
		for(int j = 0; j < optpoints.size(); j++)
		{
			glBegin(GL_LINES);
			glVertex2d(points[i].x, points[i].y);
			glVertex2d(optpoints[j].x, optpoints[j].y);
			glEnd();
		}
	}
}

bool loadMedia()
{
    //Load texture
    if( !gLoadedTexture.loadTextureFromFile( "Geometric-Modelling/cube_small.jpg" ) )
    {
        printf( "Unable to load file texture!\n" );
        return false;
    }

    return true;
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);

	GLfloat x = ( winWidth - gLoadedTexture.textureWidth() ) / 2.f;
    GLfloat y = ( winHeight - gLoadedTexture.textureHeight() ) / 2.f;
	gLoadedTexture.render( x, y );

	if(points.size() == 7)
		displayLines();
	displayPoints();
    if(optimized) {
		displayOptPoints();
		displayOptLines();
	}

	glutSwapBuffers();
}


double thesis_fn(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data){
 
    // evaluate function
    double t1 = vals_inp(0);
    double t2 = vals_inp(1);
    double t3 = vals_inp(2);
 
    double obj_val = thesisfunc(t1, t2, t3, points[0], points[4], points[5], points[6]);
 
    // update gradient
    if(grad_out){
        
        
        std::vector<double> grads = thesisfuncgradient(t1, t2, t3, points[0], points[4], points[5], points[6]);
 
        // partial derivate, t1
        (*grad_out)(0) = grads[0];
 
        // partial derivate, t2
        (*grad_out)(1) = grads[1];
        
        // partial derivate, t3
        (*grad_out)(2) = grads[2];
 
    }
 
    // return value
    return obj_val;
 
}


void optimizeForT(){
    
    std::cout << "Performing optimization..." << std::endl; 
      
    arma::vec t = arma::zeros(3,1) + 0.5;
    std::cout << "Initialized optimziation at \n" << t << std::endl;

	optim::algo_settings_t config;
	config.vals_bound = true;
	config.lower_bounds << 0 << 0 << 0;
	config.upper_bounds << 1 << 1 << 1;
      
    bool success = optim::bfgs(t, thesis_fn, nullptr, config);

    if(success){
        std::cout << "Function optimization completed successfully!" << std::endl;
    }else{
        std::cout << "Function optimization completed unsuccessfully!" << std::endl;
    }
    
    std::cout << "Solution :\n" << t << std::endl;
    
    calcOptPoints(t(0), t(1), t(2));
      
}

void keyPressed(unsigned char key, int x, int y) {
	if (key == 's') {
		optimizeForT();
	}
}

int main(int argc, char** argv) {
    

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(winWidth, winHeight);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Cube'n stuff");
	init();

	printf( "1\n" );

	//Load media
    if( !loadMedia() )
    {
        printf( "Unable to load media!\n" );
        return 2;
    }

	printf( "2\n" );

	glutDisplayFunc(display);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutKeyboardFunc(keyPressed);
	//glutSpecialFunc(keyPressed);

	update(0);
	glutMainLoop();

	return 0;
}
