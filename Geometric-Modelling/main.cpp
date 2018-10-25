#include <GL/glut.h>
#include <bevgrafmath2017.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <string>

const int winWidth = 800, winHeight = 600;
std::vector<vec2> points;

int dragged = -1;
int pointRadius = 4;
int uDist = winWidth/4.0;


void init() {

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, winWidth, 0.0, winHeight);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
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


double optFunc(double t1, double t2, double t3, double tanAlpha, double tanBeta, double tanGamma) {
	
	double sum = 0.0;

	sum += std::abs( (pow(1 - t1, 2) / pow(t1, 2)) * (pow(t2, 2) / pow(1 - t2, 2)) - tanAlpha / tanBeta );
	sum += std::abs( (pow(1 - t1, 2) / pow(t1, 2)) * (pow(t3, 2) / pow(1 - t3, 2)) - tanAlpha / tanGamma);
	sum += std::abs( (pow(1 - t2, 2) / pow(t2, 2)) * (pow(t3, 2) / pow(1 - t3, 2)) - tanBeta / tanGamma);


	return sum;

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

void display() {
	glClear(GL_COLOR_BUFFER_BIT);

	
	if(points.size() == 7)
		displayLines();
	displayPoints();

	glutSwapBuffers();
}

void keyPressed(unsigned char key, int x, int y) {
	if (key == 's') {
		//TODO
		//optimization comes here
	}
}

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(winWidth, winHeight);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Cube'n stuff");
	init();

	glutDisplayFunc(display);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	glutKeyboardFunc(keyPressed);
	//glutSpecialFunc(keyPressed);

	update(0);
	glutMainLoop();

	return 0;
}