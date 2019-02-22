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

#include "opencv2/opencv.hpp"
#include <opencv2/highgui.hpp>
std::vector<cv::Point2f> destinationPointsTmp;
const int winWidth = 800, winHeight = 600;
std::vector<vec2> points, optpoints;

// helper points for calculating "endless distant" points
// this vector will countaint 2 elements for calculating the "endless distant" points
// the first one is for Ux and Uz
// the second one is for Uy and Uz
std::vector<vec2> helperPoints;

int dragged = -1;
int fixed = -1;
int pointRadius = 4;
int uDist = winWidth/4.0;
bool optimized = false;

//File loaded texture
LTexture gLoadedTexture;
cv::Mat imageMat;
bool displayImage = true;

// default image to load
const std::string defaultImageFileName = "Geometric-Modelling/cube_small.jpg";

// variables and transformation matrices for zooming and moving
float zoom = 1.0f;
float zoomRate = 0.01f;
vec2 translateVector = vec2(0.0f);
float translateRate = 10.0f;
vec2 windowCenter = { (float) winWidth / 2.0f, (float) winHeight / 2.0f };
mat3 translateToOrigo = translate( -1 * windowCenter);
mat3 translateBack = translate(windowCenter);
mat3 scaleMatrix = scale(vec2(zoom));
mat3 translateMatrix = translate(translateVector);
mat3 M = translateBack * scaleMatrix * translateToOrigo * translateMatrix;

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

void fixPoint(std::vector<vec2> points, int sensitivity, int mouseX, int mouseY) {
    int max = points.size() > 4 ? 4 : points.size();
	for (int i = 1; i < max; i++) {
		if (dist(points[i], vec2(mouseX, winHeight - mouseY)) < sensitivity)
			fixed = i;
	}
}

double getYofLine(vec2 p1, vec2 p2, double x) {
	vec2 norm = vec2(p1.y - p2.y, -(p1.x -p2.x));

	return (norm.x * p1.x + norm.y * p1.y - norm.x * x) / norm.y;

}

double getXofLine(vec2 p1, vec2 p2, double y) {
	vec2 norm = vec2(p1.y - p2.y, -(p1.x -p2.x));

	return (norm.x * p1.x + norm.y * p1.y - norm.y * y) / norm.x;
}	

void update(int val) {
	glutPostRedisplay();
	glutTimerFunc(2, update, 0);
}

// calculate the intersection of two 2d vectors
// source: https://stackoverflow.com/a/2932601/5352042
vec2 intersect(vec2 aOrigin, vec2 aDirection, vec2 bOrigin, vec2 bDirection){

	float dx = bOrigin.x - aOrigin.x;
	float dy = bOrigin.y - aOrigin.y;

	float det = bDirection.x * aDirection.y - bDirection.y * aDirection.x;

	float u = fabs((dy * bDirection.x - dx * bDirection.y) / det);
	float v = fabs((dy * aDirection.x - dx * aDirection.y) / det);

	// calculate intersection point
	vec2 intersectionPoint = aOrigin + aDirection * u;

	return intersectionPoint;

}

// calculate "endless distant" points
void createU(){

	// base points
	vec2 o = points[0];
	vec2 x = points[1];
	vec2 y = points[2];
	vec2 z = points[3];

	// helper points
	vec2 helperPointUx = helperPoints[0];
	vec2 helperPointUy = helperPoints[1];

	// calculate direction vectors
	// o -> x
	vec2 fromOToX = (x - o);
	// z -> helperPointUx
	vec2 fromZToHelperPointUx = (helperPointUx - z);
	// o -> y
	vec2 fromOToY = (y - o);
	// z -> helperPointUy
	vec2 fromZToHelperPointUy = (helperPointUy - z);
	// x -> helperPointUx
	vec2 fromXToHelperPointUx = (helperPointUx - x);
	// o -> z
	vec2 fromOToZ = (z - o);
	// y -> helperPointUy
	vec2 fromYToHelperPointUy = (helperPointUy - y);

	// calculate Ux
	vec2 uX = intersect(o, fromOToX, z, fromZToHelperPointUx);

	// calculate Uy
	vec2 uY = intersect(o, fromOToY, z, fromZToHelperPointUy);

	// calculate Uz
	// vec2 uZ = intersect(x, fromXToHelperPointUx, o, fromOToZ);
	vec2 uZ = intersect(x, fromXToHelperPointUx, y, fromYToHelperPointUy);

	// add the calculated points to the points vector
	points.push_back(uX);
	points.push_back(uY);
	points.push_back(uZ);

}

void processMouse(int button, int action, int xMouse, int yMouse) {
	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN) {
		dragPoint(points, 10, xMouse, yMouse);

		// if we not dragging any point
		if(dragged == -1){

			// if we still picking the base points
			if(points.size() < 4){
				points.push_back(vec2(xMouse, winHeight - yMouse));
			}
			// if we still picking the helper points
			else if(helperPoints.size() != 2){

				helperPoints.push_back(vec2(xMouse, winHeight - yMouse));

				// if we have the helper points (and also the base points), lets calculate the "endless distant" points
				if(helperPoints.size() == 2){
					createU();
				}

			}

		}
			
	}
		
	if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
		dragged = -1;
    
    if(button == GLUT_RIGHT_BUTTON && action == GLUT_DOWN){
        if(fixed == -1){
            fixPoint(points, 10, xMouse, yMouse);
        }
        else
            fixed = -1;
            
    }
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
		double steepness = (points[0].x - points[dragged-3].x) / (points[0].y - points[dragged-3].y);
		
		if(std::abs(steepness) < 1){
			points[dragged].y = winHeight - yMouse;
			points[dragged].x = getXofLine(points[0], points[dragged-3], winHeight - yMouse);
		} else {
			points[dragged].x = xMouse;
			points[dragged].y = getYofLine(points[0], points[dragged-3], xMouse);
		}
	}
}

// calculate the transformed point for the display from an original point
// transformation steps: translate to origio -> scale -> translate back to center
vec2 calculateTransformedPoint(vec2 original){

	// convert original point to homogen form
	vec3 originalHomogen = ihToH(original);

	// do the transformation
	vec3 transformedHomogen = M * originalHomogen;

	// convert back the result to inhomogen form
	vec2 transformedInhomogen = hToIh(transformedHomogen);

	// return result
	return transformedInhomogen;

}

// update the global matrix used for transformation
void calculateTransformationMatrix(){

	// update scale matrix
	scaleMatrix = scale(vec2(zoom));

	// update translate matrix
	translateMatrix = translate(translateVector);

	// recalculate transformation matrix
	M = translateBack * scaleMatrix * translateToOrigo * translateMatrix;

}

void displayPoints() {

	

	// base points
	for (int i = 0; i < points.size(); i++) {
        
        if(fixed == i)
            glColor3f(0.5, 0.5, 0.5);
        else
            glColor3f(1.0, 0.0, 0.0);
        
        glPointSize(10.0);
        glBegin(GL_POINTS);

		vec2 originalPoint = points[i];
		vec2 transformedPoint = calculateTransformedPoint(originalPoint);

		glVertex2f(transformedPoint.x, transformedPoint.y);
        
        glEnd();
	}
	
	
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(10.0);
	glBegin(GL_POINTS);

	// helper points
	glColor3f(1.0, 0.0, 1.0);
	for(int i=0; i<helperPoints.size(); i++){

		vec2 originalPoint = helperPoints[i];
		vec2 transformedPoint = calculateTransformedPoint(originalPoint);

		glVertex2f(transformedPoint.x, transformedPoint.y);

	}

	glEnd();


	glColor3f(1.0, 1.0, 1.0);
	glPointSize(12.0);

	vec2 originalPoint;
	vec2 transformedPoint;

	glBegin(GL_POINTS);
	for (int i = 0; i < destinationPointsTmp.size(); i++) {

		vec2 originalPoint = vec2(destinationPointsTmp[i].x, destinationPointsTmp[i].y);
		vec2 transformedPoint = calculateTransformedPoint(originalPoint);

		glVertex2f(transformedPoint.x, transformedPoint.y);

	}
	glEnd();

}

void displayOptPoints(){

    glColor3f(0.0, 1.0, 0.0);
	glPointSize(12.0);

	vec2 originalPoint;
	vec2 transformedPoint;

	glBegin(GL_POINTS);
	for (int i = 0; i < optpoints.size(); i++) {

		vec2 originalPoint = optpoints[i];
		vec2 transformedPoint = calculateTransformedPoint(originalPoint);

		glVertex2f(transformedPoint.x, transformedPoint.y);

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

	vec2 originalPoint;
	vec2 transformedPoint;
	
	for (int i = 4; i < 7; i++) {
		glBegin(GL_LINES);

		originalPoint = points[0];
		transformedPoint = calculateTransformedPoint(originalPoint);
		glVertex2d(transformedPoint.x, transformedPoint.y);

		originalPoint = points[i];
		transformedPoint = calculateTransformedPoint(originalPoint);
		glVertex2d(transformedPoint.x, transformedPoint.y);

		glEnd();
	}

	glBegin(GL_LINE_LOOP);

	originalPoint = points[5];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	originalPoint = points[1];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	originalPoint = points[6];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	glEnd();

	glBegin(GL_LINE_LOOP);

	originalPoint = points[5];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	originalPoint = points[3];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	originalPoint = points[4];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	glEnd();

	glBegin(GL_LINE_LOOP);

	originalPoint = points[6];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	originalPoint = points[2];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	originalPoint = points[4];
	transformedPoint = calculateTransformedPoint(originalPoint);
	glVertex2d(transformedPoint.x, transformedPoint.y);

	glEnd();
	
}

void displayOptLines() {
	glColor3d(0.0, 0.0, 1.0);
	glLineWidth(3.0);

	vec2 originalPoint;
	vec2 transformedPoint;
	
	for (int i = 4; i < 7; i++) {
		glBegin(GL_LINES);

		originalPoint = points[0];
		transformedPoint = calculateTransformedPoint(originalPoint);
		glVertex2d(transformedPoint.x, transformedPoint.y);

		originalPoint = points[i];
		transformedPoint = calculateTransformedPoint(originalPoint);
		glVertex2d(transformedPoint.x, transformedPoint.y);

		glEnd();
	}

	for (int i = 4; i < 7; i++) 
	{
		for(int j = 0; j < optpoints.size(); j++)
		{
			glBegin(GL_LINES);

			originalPoint = points[i];
			transformedPoint = calculateTransformedPoint(originalPoint);
			glVertex2d(transformedPoint.x, transformedPoint.y);

			originalPoint = optpoints[j];
			transformedPoint = calculateTransformedPoint(originalPoint);
			glVertex2d(transformedPoint.x, transformedPoint.y);

			glEnd();
		}
	}
}

bool loadMedia(std::string fileName) {
    //Load texture
    if( !gLoadedTexture.loadTextureFromFile(fileName) )
    {
        printf( "Unable to load file texture!\n" );
        return false;
    }

	imageMat = cv::imread(fileName);

    return true;
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);

	if(displayImage){
		/*GLfloat x = ( winWidth - gLoadedTexture.textureWidth() ) / 2.f;
		GLfloat y = ( winHeight - gLoadedTexture.textureHeight() ) / 2.f;

		vec2 upperLeftCorner(x, y);
		vec2 upperRightCorner(x + gLoadedTexture.textureWidth(), y);
		vec2 lowerLeftCorner(x, y + gLoadedTexture.textureHeight());
		vec2 lowerRightCorner(x + gLoadedTexture.textureWidth(), y + gLoadedTexture.textureHeight());*/

		vec2 upperLeftCorner(0, 0);
		vec2 upperRightCorner(0 + gLoadedTexture.textureWidth(), 0);
		vec2 lowerLeftCorner(0, 0 + gLoadedTexture.textureHeight());
		vec2 lowerRightCorner(0 + gLoadedTexture.textureWidth(), 0 + gLoadedTexture.textureHeight());

		gLoadedTexture.render( 
			calculateTransformedPoint(upperLeftCorner), 
			calculateTransformedPoint(upperRightCorner), 
			calculateTransformedPoint(lowerLeftCorner), 
			calculateTransformedPoint(lowerRightCorner) 
		);
	}

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
    double t1, t2, t3;
    
        
    if(fixed == 1){
        t1 = (points[1].x - points[4].x) / (points[0].x - points[4].x);
        t2 = vals_inp(0);
        t3 = vals_inp(1);
        
    } else if (fixed == 2){
        t1 = vals_inp(0);
        t2 = (points[2].x - points[5].x) / (points[0].x - points[5].x);
        t3 = vals_inp(1);
        
    } else if (fixed == 3){
        t1 = vals_inp(0);
        t2 = vals_inp(1);
        t3 = (points[3].x - points[6].x) / (points[0].x - points[6].x);
        
    } else {
        t1 = vals_inp(0);
        t2 = vals_inp(1);
        t3 = vals_inp(2);
    }
    
    
 
    double obj_val = thesisfunc(t1, t2, t3, points[0], points[4], points[5], points[6]);
 
    // update gradient
    if(grad_out){
        
        std::vector<double> grads = thesisfuncgradient(t1, t2, t3, points[0], points[4], points[5], points[6]);
        
        if(fixed == 1){
            // partial derivate, t2
            (*grad_out)(0) = grads[1];
            
            // partial derivate, t3
            (*grad_out)(1) = grads[2];
        
        } else if (fixed == 2){
            // partial derivate, t1
            (*grad_out)(0) = grads[0];
            
            // partial derivate, t3
            (*grad_out)(1) = grads[2];
        
        } else if (fixed == 3){
            // partial derivate, t1
            (*grad_out)(0) = grads[0];
    
            // partial derivate, t2
            (*grad_out)(1) = grads[1];
            
        
        } else {
            // partial derivate, t1
            (*grad_out)(0) = grads[0];
    
            // partial derivate, t2
            (*grad_out)(1) = grads[1];
            
            // partial derivate, t3
            (*grad_out)(2) = grads[2];
        }
 
    }
 
    // return value
    return obj_val;
 
}


void optimizeForT(){
    
    std::cout << "Performing optimization..." << std::endl; 
    
    arma::vec t;
    if(fixed == -1)
        t = arma::zeros(3,1) + 0.5;
    else
        t = arma::zeros(2,1) + 0.5;
    
    std::cout << "Initialized optimziation at \n" << t << std::endl;

	optim::algo_settings_t config;
	config.vals_bound = true;
	config.lower_bounds << 0 << 0 << 0;
	config.upper_bounds << 1 << 1 << 1;
      
    bool success = optim::bfgs(t, thesis_fn, nullptr, config);

    
    double t1, t2, t3;
    
    if(success){
        
        if(fixed == 1){
        t1 = (points[1].x - points[4].x) / (points[0].x - points[4].x);
        t2 = t(0);
        t3 = t(1);
        
        } else if (fixed == 2){
        t1 = t(0);
        t2 = (points[2].x - points[5].x) / (points[0].x - points[5].x);
        t3 = t(1);
        
        } else if (fixed == 3){
        t1 = t(0);
        t2 = t(1);
        t3 = (points[3].x - points[6].x) / (points[0].x - points[6].x);
        
        } else {
        t1 = t(0);
        t2 = t(1);
        t3 = t(2);
        }
        
        
        std::cout.precision(15);
        double val = thesisfunc(t1, t2, t3, points[0], points[4], points[5], points[6]);
        std::cout << "Function optimization completed successfully!" << std::endl;
        std::cout << "Value of Szabo-theoriem: " << std::fixed << val << std::endl;
    }else{
        std::cout << "Function optimization completed unsuccessfully!" << std::endl;
    }
    
    std::cout << "Solution :\n" << t << std::endl;
    
    calcOptPoints(t1, t2, t3);
      
}

void scaleDownOptimizedPoints() {
	double dist1 = dist(points[1], optpoints[0]);
	double distOrigOpt1 = dist(points[0], optpoints[0]);

	double dist2 = dist(points[2], optpoints[1]);
	double distOrigOpt2 = dist(points[0], optpoints[1]);

	double dist3 = dist(points[3], optpoints[2]);
	double distOrigOpt3 = dist(points[0], optpoints[2]);

	double scaleRatio;
	bool scaleDown;

	if(dist1 <= dist2 && dist1 <= dist3) {
		scaleRatio = dist1 / distOrigOpt1;

		if(dist(points[0], points[1]) <= distOrigOpt1) {
			scaleDown = true;
		} else {
			scaleDown = false;
		}
	} else if(dist2 <= dist1 && dist2 <= dist3) {
		scaleRatio = dist2 / distOrigOpt2;
		
		if(dist(points[0], points[2]) <= distOrigOpt2) {
			scaleDown = true;
		} else {
			scaleDown = false;
		}
	} else if(dist3 <= dist1 && dist3 <= dist2) {
		scaleRatio = dist3 / distOrigOpt3;
		
		if(dist(points[0], points[3]) <= distOrigOpt3) {
			scaleDown = true;
		} else {
			scaleDown = false;
		}
	}

	if(scaleDown) {
		optpoints[0] = points[0] + normalize(optpoints[0] - points[0]) * (distOrigOpt1 - scaleRatio * distOrigOpt1);
		optpoints[1] = points[0] + normalize(optpoints[1] - points[0]) * (distOrigOpt2 - scaleRatio * distOrigOpt2);
		optpoints[2] = points[0] + normalize(optpoints[2] - points[0]) * (distOrigOpt3 - scaleRatio * distOrigOpt3);
	} else {
		optpoints[0] = points[0] + normalize(optpoints[0] - points[0]) * (distOrigOpt1 + scaleRatio * distOrigOpt1);
		optpoints[1] = points[0] + normalize(optpoints[1] - points[0]) * (distOrigOpt2 + scaleRatio * distOrigOpt2);
		optpoints[2] = points[0] + normalize(optpoints[2] - points[0]) * (distOrigOpt3 + scaleRatio * distOrigOpt3);
	}
}

//TODO im sure these are transformed points or stuff like that, so it sill sucks, or if not transformed, it needs to be scaled on to the image.
void transformImage() {

	std::vector<vec2> originalIntersections;

	vec2 fromXToUz = points[6] - points[1];
	vec2 fromZToUx = points[4] - points[3];
	originalIntersections.push_back(intersect(points[1], fromXToUz, points[3], fromZToUx));

	vec2 fromYToUz = points[6] - points[2];
	vec2 fromZToUy = points[5] - points[3];
	originalIntersections.push_back(intersect(points[2], fromYToUz, points[3], fromZToUy));

	vec2 fromXToUy = points[5] - points[1];
	vec2 fromYToUx = points[4] - points[2];
	originalIntersections.push_back(intersect(points[1], fromXToUy, points[2], fromYToUx));

	std::vector<cv::Point2f> sourcePoints;
	sourcePoints.push_back(cv::Point2f(points[0].x, imageMat.rows - points[0].y));
    sourcePoints.push_back(cv::Point2f(points[1].x, imageMat.rows - points[1].y));
    sourcePoints.push_back(cv::Point2f(originalIntersections[0].x, imageMat.rows - originalIntersections[0].y));
    sourcePoints.push_back(cv::Point2f(points[3].x, imageMat.rows - points[3].y));
	sourcePoints.push_back(cv::Point2f(originalIntersections[1].x, imageMat.rows - originalIntersections[1].y));
	sourcePoints.push_back(cv::Point2f(points[2].x, imageMat.rows - points[2].y));
	sourcePoints.push_back(cv::Point2f(originalIntersections[2].x, imageMat.rows - originalIntersections[2].y));
	
	/*std::cout << sourcePoints[0].x << " " << sourcePoints[0].y << std::endl;
	std::cout << sourcePoints[1].x << " " << sourcePoints[1].y << std::endl;
	std::cout << sourcePoints[2].x << " " << sourcePoints[2].y << std::endl;
	std::cout << sourcePoints[3].x << " " << sourcePoints[3].y << std::endl;
	std::cout << sourcePoints[4].x << " " << sourcePoints[4].y << std::endl;
	std::cout << sourcePoints[5].x << " " << sourcePoints[5].y << std::endl;*/

	std::vector<vec2> optimizedIntersections;

	vec2 fromOptimozedXToUz = points[6] - optpoints[0];
	vec2 fromOptimizedZToUx = points[4] - optpoints[2];
	optimizedIntersections.push_back(intersect(optpoints[0], fromOptimozedXToUz, optpoints[2], fromOptimizedZToUx));

	vec2 fromOptimizedYToUz = points[6] - optpoints[1];
	vec2 fromOptimizedZToUy = points[5] - optpoints[2];
	optimizedIntersections.push_back(intersect(optpoints[1], fromOptimizedYToUz, optpoints[2], fromOptimizedZToUy));

	vec2 fromOptimizedXToUy = points[5] - optpoints[0];
	vec2 fromOptimizedYToUx = points[4] - optpoints[1];
	optimizedIntersections.push_back(intersect(optpoints[0], fromOptimizedXToUy, optpoints[1], fromOptimizedYToUx));

	std::vector<cv::Point2f> destinationPoints;
	destinationPoints.push_back(cv::Point2f(points[0].x, imageMat.rows - points[0].y));
    destinationPoints.push_back(cv::Point2f(optpoints[0].x, imageMat.rows - optpoints[0].y));
    destinationPoints.push_back(cv::Point2f(optimizedIntersections[0].x, imageMat.rows - optimizedIntersections[0].y));
    destinationPoints.push_back(cv::Point2f(optpoints[2].x, imageMat.rows - optpoints[2].y));
    destinationPoints.push_back(cv::Point2f(optimizedIntersections[1].x, imageMat.rows - optimizedIntersections[1].y));
    destinationPoints.push_back(cv::Point2f(optpoints[1].x, imageMat.rows - optpoints[1].y));
    destinationPoints.push_back(cv::Point2f(optimizedIntersections[2].x, imageMat.rows - optimizedIntersections[2].y));


	destinationPointsTmp.push_back(cv::Point2f(points[0].x, points[0].y));
    destinationPointsTmp.push_back(cv::Point2f(optpoints[0].x, optpoints[0].y));
    destinationPointsTmp.push_back(cv::Point2f(optimizedIntersections[0].x, optimizedIntersections[0].y));
    destinationPointsTmp.push_back(cv::Point2f(optpoints[2].x, optpoints[2].y));
    destinationPointsTmp.push_back(cv::Point2f(optimizedIntersections[1].x, optimizedIntersections[1].y));
    destinationPointsTmp.push_back(cv::Point2f(optpoints[1].x, optpoints[1].y));
    destinationPointsTmp.push_back(cv::Point2f(optimizedIntersections[2].x, optimizedIntersections[2].y));



	cv::Mat homography = findHomography(sourcePoints, destinationPoints, 0, 1);

	cv::Mat transformedImageMat;

	warpPerspective(imageMat, transformedImageMat, homography, imageMat.size());

	/*GLuint* imageArray = (GLuint*)calloc(transformedImageMat.rows * transformedImageMat.cols, sizeof(GLuint));
	std::cout << "asd" << std::endl;
	for(int i = 0; i < transformedImageMat.rows; ++i) {
		for(int j = 0; j < transformedImageMat.cols; ++j) {
			imageArray[i * transformedImageMat.cols + j] = (GLuint)transformedImageMat.ptr<int>(i)[j];
			std::cout << i << " " << j << std::endl;
		}
	}
	std::cout << "asd" << std::endl;*/

	imwrite("newCube.jpg", transformedImageMat);

	gLoadedTexture.freeTexture();

	if( !gLoadedTexture.loadTextureFromFile("newCube.jpg") )
    {
        printf( "Unable to load the transformed file texture!\n" );
        return;
    }

	//gLoadedTexture.loadTextureFromPixels32(imageArray, transformedImageMat.cols, transformedImageMat.rows);

	//imshow("Source Image", imageMat);
    //imshow("Warped Source Image", transformedImageMat);
 
    //cv::waitKey(0);
}

void keyPressed(unsigned char key, int x, int y) {
	if (key == 's') {
		optimizeForT();
	} else if (key == 'd') {
		scaleDownOptimizedPoints();
	} else if (key == 'f') {
		transformImage();
	} else if (key == 'q') {
		zoom -= zoomRate;
		calculateTransformationMatrix();
	} else if (key == 'w') {
		zoom += zoomRate;
		calculateTransformationMatrix();
	} else if (key == 'r') {

		// reset zoom
		zoom = 1.0f;

		// reset translate
		translateVector = vec2(0.0f);

		calculateTransformationMatrix();

	} else if (key == 'x') {
		// turn on/off image rendering
		displayImage = !displayImage;
	} else if(key == 'i') {
		translateVector.y -= translateRate;
		calculateTransformationMatrix();
	} else if(key == 'j') {
		translateVector.x += translateRate;
		calculateTransformationMatrix();
	} else if(key == 'k') {
		translateVector.y += translateRate;
		calculateTransformationMatrix();
	} else if(key == 'l') {
		translateVector.x -= translateRate;
		calculateTransformationMatrix();
	}
}

int main(int argc, char** argv) {
    

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(winWidth, winHeight);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Cube'n stuff");
	init();

	std::string fileName = defaultImageFileName;
	if(argc == 2){
		fileName = argv[1];
	}

	//Load media
	if( !loadMedia(fileName) )
	{
		printf( "Unable to load media!\n" );
		return 2;
	}

	glutDisplayFunc(display);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutKeyboardFunc(keyPressed);
	//glutSpecialFunc(keyPressed);

	update(0);
	glutMainLoop();

	return 0;
}
