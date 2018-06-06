//
//  main.cpp
//  CG lab2
//
//  Created by Janice Wang on 3/1/18.
//  Copyright © 2018 Janice Wang. All rights reserved.
//

#include <OpenGL/gl.h>
#include <GLUT/GLUT.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
using namespace std;

class vector3d
{
    
public:
    float x;
    float y;
    float z;
    vector3d(float posX, float posY, float posZ){
        x=posX;
        y=posY;
        z=posZ;
    }
    vector3d(){
        x=0;
        y=0;
        z=0;
    }
    void setCoord(float newX, float newY, float newZ){
        x=newX;
        y=newY;
        z=newZ;
    }
    
    vector3d operator*(vector3d v){
        vector3d res;
        float tempX = y*v.z - z*v.y;
        float tempY = z*v.x - x*v.z;
        float tempZ = x*v.y - y*v.x;
        res.setCoord(tempX, tempY, tempZ);
        return res;
    }
    
    vector3d operator-(vector3d v){
        vector3d res;
        float tempX = x - v.x;
        float tempY = y - v.y;
        float tempZ = z - v.z;
        res.setCoord(tempX, tempY, tempZ);
        return res;
    }
    
    vector3d operator/(float n){
        vector3d res;
        float tempX = x / n;
        float tempY = y / n;
        float tempZ = z / n;
        res.setCoord(tempX, tempY, tempZ);
        return res;
    }
    
    float length(){
        float res;
        float temp = pow(x,2) + pow(y,2) + pow(z,2);
        res=sqrt(temp);
        return res;
    }
    
};

//view transaction. C is camera position, P is view point position, Vp is V prime

vector3d viewTrans(vector3d C, vector3d P, vector3d Vp,vector3d vetex){
    vector3d vp;
    vector3d tempN = P - C;
    vector3d N = tempN / tempN.length();
    vector3d tempU = N * Vp;
    vector3d U = tempU / tempU.length();
    vector3d V = U * N;
    
    float Mview[4][4]={{U.x,U.y,U.z,-(U.x*C.x+U.y*C.y+U.z*C.z)},
        {V.x,V.y,V.z,-(V.x*C.x+V.y*C.y+V.z*C.z)},
        {N.x,N.y,N.z,-(N.x*C.x+N.y*C.y+N.z*C.z)},
        {0,0,0,1}};
    float world[4]={vetex.x,vetex.y,vetex.z,1};
    float temp[4];
    float Vcoord;
    for(int j=0;j<4;j++){
        Vcoord=0;
        for(int i=0;i<4;i++){
            Vcoord+=Mview[j][i]*world[i];
        }
        temp[j]=Vcoord;
    }
    vp.setCoord(temp[0], temp[1], temp[2]);
    return vp;
}

//create Mpers
vector3d persTrans(float d, float h,float f,vector3d vertex){
    vector3d res;
    float px=d * vertex.x / h * vertex.z;
    float py=d * vertex.y / h * vertex.z;
    float pz=f * (1-d/vertex.z) / (f-d);
    res.setCoord(px, py, pz);
    return res;
}

//scal object

vector3d scal(vector3d v){
    vector3d res;
    int s=2000;
    res.x = v.x*s;
    res.y = v.y*s;
    res.z = v.z;
    return res;
}

//Vnum is vertices number, Pnum is the number of polygon.

int Vnum,Pnum;
vector3d *v;
int **p;//polygon

//Scan vertices and polygon from .d file. For z-buffer algorithm,

void read()
{
    
    std::ifstream in("house.d");
    string data;
    in >> data >> Vnum >> Pnum;
    v = new vector3d[2*Vnum];
    
    for(int i = 0; i < Vnum; i++){
        in>> v[i].x >> v[i].y >>v[i].z;
        int Vn = i + Vnum;
        v[Vn].x = v[i].x + 20;
        v[Vn].y = v[i].y;
        v[Vn].z = v[i].z + 10;
    }
    
    vector3d C,P,Vp;
    C.setCoord(30,50,15);
    P.setCoord(20,20,50);
    Vp.setCoord(0,1,0);
    
    float d,h,f;
    d = 0.2;
    h = 500;
    f = 6;
    
    for(int i = 0; i < 2 * Vnum; i++){
        v[i] = viewTrans(C, P, Vp, v[i]);
        v[i] = persTrans(d, h, f, v[i]);
        v[i] = scal(v[i]);
    }
    
    p = new int*[2*Pnum];
    for(int i = 0; i < Pnum; i++){
        int e;//edges
        in >> e;
        int newE = i + Pnum;
        p[i] = new int[e+1];
        p[i][0] = e;
        
        for(int j = 1;j < e+1; j++){
            in >> p[i][j];
        }
        
        p[newE] = new int[e+1];
        p[newE][0] = e;
        
        for(int j = 1; j < e+1;j++){
            p[newE][j] = p[i][j] + Vnum;
        }
    }
}



//Each edge is made of two vertexs.
//The class contains 4 functions of calculating the minimal of x,y ,max of y and slope of edge

class Edge
{
    
public:
    int a;
    int b;
    Edge(){
        a=-1;
        b=-1;
    }
    void setEdge(int s,int q){
        a=s;
        b=q;
    }
    
    float Xmin(){
        return (v[a].y<v[b].y) ? v[a].x : v[b].x;
    }
    
    float Ymin(){
        return (v[a].y>v[b].y) ? v[b].y : v[a].y;
    }
    
    float Ymax(){
        return (v[a].y>v[b].y) ? v[a].y : v[b].y;
    }
    
    float Slop(){
        return (v[b].x-v[a].x)/(v[b].y-v[a].y);
    }
};

//Structure used to store data of pixel. It includes x, y coordinate, depth and color of pixel.
struct pixel {
    int x;
    int y;
    float depth;
    int i;
};

//find minimal scan line

void findSL(int * SL){
    for(int i = 0; i < 2*Pnum; i++){
        int Enum=p[i][0];
        
        vector3d v1=v[p[i][1]-1]-v[p[i][2]-1];
        vector3d v2=v[p[i][3]-1]-v[p[i][2]-1];
        vector3d v3=v1*v2;
        vector3d nv = v3/v3.length();
        
        if(nv.z <= 0){
            for(int j = 1; j < Enum + 1;j++){
                if(v[p[i][j]-1].y<*SL){
                    *SL=(int)v[p[i][j]-1].y;
                }
                
            }
        }
    }
}

//store information of edge table and active edge table.
//Or it stores information of edges. Ymin and Ymax In addition to z and y coordinates of two vertexs.
struct element{
    int Ymin;
    int Ymax;
    float X;
    float m;
    float z1;
    float z2;
    float y1;
    float y2;
};

//The class used to store edge table and active edge table of each polygon.
//Constructor of class converts each polygon into edge table.
//Then fill in all pixels in scan line. Furthermore, color of polygon is also stored.

class Polygon
{
    
public:
    vector<element> ET;
    vector<element> ActiveET;
    int miny=INT_MAX;
    int color;
    Polygon(int z, int color){
        this->color = color;
        vector<Edge> E;
        int Enum = p[z][0];
        for(int i = 1; i < Enum+1; i++){
            Edge newE;
            if(i == Enum){
                newE.setEdge(p[z][i]-1,p[z][1]-1);
            }else{
                newE.setEdge(p[z][i]-1,p[z][i+1]-1);
            }
            if(v[p[z][i]-1].y<miny){
                miny=(int)v[p[z][i]-1].y;
            }
            E.push_back(newE);
        }
        
        for(int i = 0; i < E.size(); i++){
            Edge temp=E[i];
            float z1,z2,y1,y2;
            z1=v[temp.a].z;
            z2=v[temp.b].z;
            y1=v[temp.a].y;
            y2=v[temp.b].y;
            
            if(v[temp.a].y==v[temp.b].y){
                float maxx,minx;
                
                if(v[temp.a].x>v[temp.b].x){
                    maxx=v[temp.a].x;
                    minx=v[temp.b].x;
                }else{
                    maxx=v[temp.b].x;
                    minx=v[temp.a].x;
                }
                
                element e1={(int)temp.Ymin(),(int)temp.Ymax(),minx,0,z1,z2,y1,y2};
                element e2={(int)temp.Ymin(),(int)temp.Ymax(),maxx,0,z1,z2,y1,y2};
                ET.push_back(e1);
                ET.push_back(e2);
            }else{
                element Newe={(int)temp.Ymin(),(int)temp.Ymax(),temp.Xmin(),temp.Slop(),z1,z2,y1,y2};
                ET.push_back(Newe);
            }
        }
    }
    
    static bool compare(element a, element b){
        return a.X < b.X;
    }
    bool scanCon(int sl,map<string,pixel>*zbuffer){
        if(sl<miny){
            return true;
        }
        if(!ET.empty()||!ActiveET.empty()){
            vector<element>::iterator it;
            for(it = ET.begin();it!= ET.end();){
                if((*it).Ymin==sl){
                    ActiveET.push_back((*it));
                    it=ET.erase(it);
                }else{
                    it++;
                }
            }
            sort(ActiveET.begin(), ActiveET.end(), compare);
            
            for(int i=0; i < ActiveET.size()-1; i+=2){
                
                float za,zb;
                float z1 = ActiveET[i].z1;
                float z2 = ActiveET[i].z2;
                float y1 = ActiveET[i].y1;
                float y2 = ActiveET[i].y2;
                
                za = z1-(z1-z2)*(y1-sl)/(y1-y2);
                zb = z1-(z1-z2)*(y1-sl)/(y1-y2);
                
                z1 = ActiveET[i+1].z1;
                z2 = ActiveET[i+1].z2;
                y1 = ActiveET[i+1].y1;
                y2 = ActiveET[i+1].y2;
                
                
                int upperXi=(int)ActiveET[i].X;
                int lowerXi=(int)(ActiveET[i+1].X)+1;
                
                if(lowerXi-1 == upperXi+1){
                    lowerXi = upperXi;
                }
                for(int j=upperXi;j<=lowerXi;j++){
                    string key = to_string(sl) + "," + to_string(j);
                    float z = zb - (zb-za)*(lowerXi-j)/(lowerXi-upperXi);
                    pixel newpixel = {j,sl,z,color};
                    
                    if((*zbuffer).empty()){
                        (*zbuffer).insert(pair<string, pixel>(key,newpixel));
                    }else{
                        map<string,pixel>::iterator zit;
                        zit=(*zbuffer).find(key);
                        if(zit==(*zbuffer).end()){
                            (*zbuffer).insert(pair<string, pixel>(key,newpixel));
                        }else{
                            if(z>zit->second.depth){
                                zit->second=newpixel;
                            }
                        }
                    }
                }
            }
            vector<element>::iterator newit;
            for(newit = ActiveET.begin(); newit != ActiveET.end();){
                if((*newit).Ymax == sl){
                    newit = ActiveET.erase(newit);
                }else{
                    (*newit).X+=(*newit).m;
                    newit++;
                }
            }
            sort(ActiveET.begin(), ActiveET.end(), compare);
            return true;
        }else{
            return false;
        }
    }
};

//draw objects in view port.
//First, store data of each polygon in 'Polygon' and distridute color to each 'Polygon'.
//Then scan all polygons, get all pixel in scan line and store pixels in 'zbuffer'.
//Finally, draw pixels and clear 'zbuffer' for next scan line 'scanLine++'.
void drawPixel(void){
    read();
    int color[7][3]={{1,0,0},{0,1,0},{0,0,1},{1,0,1},{1,1,0},{0,1,1},{1,1,1}};
    int colorIndex=0;
    glOrtho(-2000,2000,-2000,2000,-2000,2000);
    glClear(GL_COLOR_BUFFER_BIT);
    GLfloat pointSize = 1.0f;
    glPointSize(pointSize);
    map<string,pixel> zbuffer;
    int scanLine;
    findSL(&scanLine);
    vector<Polygon> Ps;
    for(int i=0;i<2*Pnum;i++){
        vector3d v1=v[p[i][1]-1]-v[p[i][2]-1];
        vector3d v2=v[p[i][3]-1]-v[p[i][2]-1];
        vector3d v3=v1*v2;
        vector3d nv = v3/v3.length();
        if(nv.z<=0){
            Polygon newp(i,colorIndex);
            Ps.push_back(newp);
            if(colorIndex==6){
                colorIndex=0;
            }else{
                colorIndex++;
            }
        }
    }
    while(!Ps.empty()){
        vector<Polygon>::iterator pit;
        for(pit=Ps.begin();pit!=Ps.end();){
            if((*pit).scanCon(scanLine,&zbuffer)){
                pit++;
            }else{
                Ps.erase(pit);
            }
        }
        
        glBegin(GL_POINTS);
        map<string,pixel>::iterator zit;
        for(zit=zbuffer.begin();zit!=zbuffer.end();){
            pixel ptemp=zit->second;
            glColor3f(color[ptemp.i][0], color[ptemp.i][1], color[ptemp.i][2]);
            glVertex2d(ptemp.x, ptemp.y);
            zit++;
        }
        glEnd();
        zbuffer.clear();
        scanLine++;
    }
    glFlush();
}

int main(int argc, char** argv){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(640,480);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("CG lab2");
    glutDisplayFunc(drawPixel);
    glutMainLoop();
    return 0;
}



