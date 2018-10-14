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


void init() {

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	gluOrtho2D(0.0, winWidth, 0.0, winHeight);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
}


void dragPoint(std::vector<vec2> points, int sensitivity, int mouseX, int mouseY) {
	for (int i = 0; i < points.size(); i++) {
		if (dist(points[i], vec2(mouseX, winHeight - mouseY)) < sensitivity)
			dragged = i;
	}
}

void update(int val) {
	glutPostRedisplay();
	glutTimerFunc(2, update, 0);
}

void processMouse(int button, int action, int xMouse, int yMouse) {
	if (button == GLUT_LEFT_BUTTON && action == GLUT_DOWN) {
		dragPoint(points, 10, xMouse, yMouse);
		if (dragged == -1)
			points.push_back(vec2(xMouse, winHeight - yMouse));
	}
		
	if (button == GLUT_LEFT_BUTTON && action == GLUT_UP)
		dragged = -1;
}

void processMouseActiveMotion(int xMouse, int yMouse) {
	if (dragged >= 0) {
		points[dragged].x = xMouse;
		points[dragged].y = winHeight - yMouse;
	}
}

void displayPoints() {

	glColor3f(0.0, 0.0, 0.0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	for (int i = 0; i < points.size(); i++) {
		glVertex2f(points[i].x, points[i].y);
	}
	glEnd();

}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);

	displayPoints();

	glutSwapBuffers();
}

void keyPressed(int key, int x, int y) {
	
}

int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(winWidth, winHeight);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Curves'n stuff");
	init();

	glutDisplayFunc(display);
	glutMouseFunc(processMouse);
	glutMotionFunc(processMouseActiveMotion);
	//glutKeyboardFunc(keyPressed);
	glutSpecialFunc(keyPressed);

	update(0);
	glutMainLoop();

	return 0;
}