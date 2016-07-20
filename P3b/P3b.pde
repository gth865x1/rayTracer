//Catherine Chapman
//3451

// This is the starter code for the CS 3451 Ray Tracing project.
// The most important part of this code is the interpreter, which will
// help you parse the scene description (.cli) files.

// A global variable for holding current active file name.
// By default, the program reads in i0.cli, which draws a rectangle.
String gCurrentFile = new String("i0.cli");
float fov, Car, Cag, Cab, Cdr, Cdg, Cdb, Csr, Csg, Csb, sphereR, sphereX, sphereY, sphereZ, 
  lightR, lightG, lightB, lightX, lightY, lightZ, r, g, b, P, Krefl, vertX, vertY, vertZ;
int depth;
ArrayList<Objects> objects = new ArrayList();
ArrayList<Light> lights = new ArrayList();
ArrayList<Float> vertices = new ArrayList();

Sphere newSphere;
Triangle newTriangle;
Light newLight;
boolean isHit;


void setup() {
  size(500, 500);  
  noStroke();
  colorMode(RGB, 1.0);
  background(0, 0, 0);
  interpreter();
}

void reset_scene() {
  //reset global scene variables
  background(0, 0, 0);
  objects.clear();
  lights.clear();  
  vertices.clear();
  depth = 0;
}

void keyPressed() {
  reset_scene();
  switch(key) {
  case '1':  
    gCurrentFile = new String("i1.cli"); 
    interpreter(); 
    break;
  case '2':  
    gCurrentFile = new String("i2.cli"); 
    interpreter(); 
    break;
  case '3':  
    gCurrentFile = new String("i3.cli"); 
    interpreter(); 
    break;
  case '4':  
    gCurrentFile = new String("i4.cli"); 
    interpreter(); 
    break;
  case '5':  
    gCurrentFile = new String("i5.cli"); 
    interpreter(); 
    break;
  case '6':  
    gCurrentFile = new String("i6.cli"); 
    interpreter(); 
    break;
  case '7':  
    gCurrentFile = new String("i7.cli"); 
    interpreter(); 
    break;
  case '8':  
    gCurrentFile = new String("i8.cli"); 
    interpreter(); 
    break;
  case '9':  
    gCurrentFile = new String("i9.cli"); 
    interpreter(); 
    break;
  case '0':  
    gCurrentFile = new String("i10.cli"); 
    interpreter(); 
    break;
  }
}

float get_float(String str) { 
  return float(str);
}

// this routine helps parse the text in a scene description file
void interpreter() {
  println("Parsing '" + gCurrentFile + "'");
  String str[] = loadStrings(gCurrentFile);
  if (str == null) println("Error! Failed to read the file.");
  for (int i=0; i<str.length; i++) {

    String[] token = splitTokens(str[i], " "); // Get a line and parse tokens.
    if (token.length == 0) continue; // Skip blank line.

    if (token[0].equals("fov")) {
      fov = get_float(token[1]);
    } else if (token[0].equals("background")) {
      r =get_float(token[1]);
      g =get_float(token[2]);
      b =get_float(token[3]);
    } else if (token[0].equals("light")) {
      lightX =get_float(token[1]);
      lightY =get_float(token[2]);
      lightZ =get_float(token[3]);
      lightR =get_float(token[4]);
      lightG =get_float(token[5]);
      lightB =get_float(token[6]);
      newLight = new Light(lightR, lightG, lightB, lightX, lightY, lightZ);
      lights.add(newLight);
      // println("Lights: "+ lights.size());
    } else if (token[0].equals("surface")) {
      Cdr =get_float(token[1]);
      Cdg =get_float(token[2]);
      Cdb =get_float(token[3]);
      Car =get_float(token[4]);
      Cag =get_float(token[5]);
      Cab =get_float(token[6]);
      Csr =get_float(token[7]);
      Csg =get_float(token[8]);
      Csb =get_float(token[9]);
      P =get_float(token[10]);
      Krefl =get_float(token[11]);
    } else if (token[0].equals("sphere")) {
      sphereR =get_float(token[1]);
      sphereX =get_float(token[2]);
      sphereY =get_float(token[3]);
      sphereZ =get_float(token[4]);
      newSphere = new Sphere(Cdr, Cdg, Cdb, Csr, Csg, Csb, P, Krefl, Car, Cag, Cab);
      objects.add(newSphere);
    } else if (token[0].equals("begin")) {
      newTriangle = new Triangle(Cdr, Cdg, Cdb, Csr, Csg, Csb, P, Krefl, Car, Cag, Cab);
    } else if (token[0].equals("vertex")) {
      vertX =get_float(token[1]);
      vertY =get_float(token[2]);
      vertZ =get_float(token[3]);
      vertices.add(vertX);
      vertices.add(vertY);
      vertices.add(vertZ);
    } else if (token[0].equals("end")) {

      newTriangle.p1.x = vertices.get(0);
      newTriangle.p1.y = vertices.get(1);
      newTriangle.p1.z = vertices.get(2);
      newTriangle.p2.x = vertices.get(3);
      newTriangle.p2.y = vertices.get(4);
      newTriangle.p2.z = vertices.get(5);
      newTriangle.p3.x = vertices.get(6);
      newTriangle.p3.y = vertices.get(7);
      newTriangle.p3.z = vertices.get(8);
      objects.add(newTriangle);
    } else if (token[0].equals("rect")) {   // this command demonstrates how the parser works
      float x =get_float(token[1]);       // and is not really part of the ray tracer
      float y =get_float(token[2]);
      float w =get_float(token[3]);
      float h =get_float(token[4]);
      fill (255, 255, 255);  // make the fill color white
      rect (x, y, w, h);     // draw a rectangle on the screne
    } else if (token[0].equals("write")) {
      draw_scene();   // this is where you actually perform the ray tracing
      println("Saving image to '" + token[1]+"'");
      save(token[1]); // this saves your ray traced scene to a PNG file
    }
  }
}

// This is where you should put your code for creating
// eye rays and tracing them.
void draw_scene() {
  depth=0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      //if (x!=y) continue;
      // create and cast an eye ray
      float k = tan(radians(fov)/2);
      float xPrime = (x-(width/2))*((2*k)/width);
      float yPrime = -((y-(height/2))*((2*k)/height));
      float zPrime = -1;
      PVector colorHere = new PVector();
      Ray ray = new Ray(0, 0, 0, xPrime, yPrime, zPrime);
      colorHere = colorRay(ray);
      //set the pixel color
      fill (colorHere.x, colorHere.y, colorHere.z);// you should put the correct pixel color here
      rect (x, y, 1, 1);  // make a tiny rectangle to fill the pixel
    }
  }
}

void draw() {
  // nothing should be placed here for this project
}

PVector colorRay(Ray ray) {

  PVector rayColor = new PVector();  
  Hit firstHit = null;
  float tMin = Float.MAX_VALUE;
  int hitObject = -1;
  for (int i=0; i<objects.size(); i++) {
    Hit boom = raySceneInteraction(ray, i);
    if ((boom.time > 0) && (boom.time < tMin)) {
      tMin = boom.time;
      hitObject = i;
      firstHit = new Hit(tMin, ray, hitObject, boom.normal);
    }
  }
  if (firstHit == null) {
    rayColor.x = r;
    rayColor.y = g;
    rayColor.z = b;
    return rayColor;
  } else {
    rayColor = shadeHit(firstHit, hitObject, depth);
    return rayColor;
  }
}

Hit raySceneInteraction(Ray ray, int i) {

  float t = Float.MAX_VALUE;
  float t1 = 0;
  float t2 = 0;

  float dx = ray.direction.x;
  float dy = ray.direction.y;
  float dz = ray.direction.z;

  float origX = ray.origin.x;
  float origY = ray.origin.y;
  float origZ = ray.origin.z;

  int hitObNum = 0;

  PVector hitNormal = new PVector();
  float hitX, hitY, hitZ;
  PVector hitPoint = new PVector();
  PVector center;

  isHit = true;
  Hit boom;

  if (objects.get(i).objectID ==1) {
    Sphere sphere = (Sphere)objects.get(i);
    float centerX = sphere.sphereCenter.x;
    float centerY = sphere.sphereCenter.y;
    float centerZ = sphere.sphereCenter.z;
    float sphereRad = sphere.radius;

    center = new PVector(centerX, centerY, centerZ);

    float a = sq(dx) + sq(dy) + sq(dz);
    float b = 2*(((origX-centerX)*dx)+((origY-centerY)*dy)+((origZ-centerZ)*dz));

    float c = sq(origX-centerX)+sq(origY-centerY)+sq(origZ-centerZ) - sq(sphereRad);

    float underTheRad = sq(b) - (4*a*c);


    if (underTheRad >= 0) {
      t1 = (-b + sqrt(underTheRad))/(2*a);
      t2 = (-b - sqrt(underTheRad))/(2*a);

      if (t1 < 0 && t2>0) {
        t = t2;
      }
      if (t2 < 0 && t1>0) {
        t= t1;
      }
      if (t1 < t2) {
        if (t1 < t) {
          t = t1;
          isHit = true;
          hitObNum = i;
        }
      } else if (t1 > t2) {
        if (t2 < t) {
          t = t2;
          isHit = true;
          hitObNum = i;
        }
      } else if (t1 == t2) {
        if (t1 < t) {
          t = t1;
          isHit = true;
          hitObNum = i;
        }
      }
      hitX = origX+(t*dx);
      hitY = origY+(t*dy);
      hitZ = origZ+(t*dz);
      hitPoint.x = hitX;
      hitPoint.y = hitY;
      hitPoint.z = hitZ;

      PVector.sub(hitPoint, center, hitNormal);
      hitNormal.normalize();
    } else {
      t=-1;
      isHit = false;
    }
  } else if (objects.get(i).objectID ==2) {
    //triangle/polygon intersection
    PVector e1, e2, norm, eCross;
    Triangle tri = (Triangle)objects.get(i);
    e1 = PVector.sub(tri.p2, tri.p1);
    e2 = PVector.sub(tri.p3, tri.p1);

    eCross = new PVector();
    PVector.cross(e1, e2, eCross);

    norm = eCross.normalize(null);

    float a = norm.x;
    float b = norm.y;
    float c = norm.z;

    float d = -((a*tri.p1.x)+(b*tri.p1.y)+(c*tri.p1.z));

    if (norm.dot(ray.direction) == 0) {
      t=-1;
    } else {
      float tempT = -((a*origX)+(b*origY)+(c*origZ)+d)/((a*dx)+(b*dy)+(c*dz));

      float x = origX+tempT*dx;
      float y = origY+tempT*dx;
      float z = origZ+tempT*dz;

      PVector hit = new PVector(x, y, z);

      PVector v1 = PVector.sub(tri.p2, tri.p1); 
      PVector v2 = PVector.sub(tri.p3, tri.p2);
      PVector v3 = PVector.sub(tri.p1, tri.p3);
      v1.normalize();
      v2.normalize();
      v3.normalize();

      PVector hitVector1 = PVector.sub(hit, tri.p1);
      PVector hitVector2 = PVector.sub(hit, tri.p2);
      PVector hitVector3 = PVector.sub(hit, tri.p3);

      float v1Hit = v1.cross(hitVector1).dot(norm);
      float v2Hit = v2.cross(hitVector2).dot(norm);
      float v3Hit = v3.cross(hitVector3).dot(norm);

      if (v1Hit >0 && v2Hit >0 && v3Hit > 0) {
        if (tempT < t) {
          t = tempT;
          hitObNum = i;
        }
      }
    }
    hitNormal = norm;
    //    println("this is t: " + t);
  }

  boom = new Hit(t, ray, hitObNum, hitNormal);

  return boom;
}

PVector shadeHit(Hit hit, int obj, int depth) {
  if (depth > 10) {
    return hit.hitColor;
  } else {
    PVector hitPoint = hit.hitPoint;

    PVector eyeRay = hit.eyeRay.direction;
    PVector.mult(eyeRay, -1, eyeRay);

    if (objects.get(obj).objectID == 2) {
      PVector.mult(hit.normal, -1, hit.normal);
    }


    float visibleLi = 0;
    PVector lightPos;
    PVector lightVec = new PVector();

    float lightColorX = 0;
    float lightColorY = 0;
    float lightColorZ = 0;

    float phongHL;

    PVector reflection;
    float normdotEye = hit.normal.dot(eyeRay);
    float normMath = 2*normdotEye;
    PVector normThird = PVector.mult(hit.normal, normMath);
    reflection = PVector.sub(normThird, eyeRay);
    reflection.normalize();
    Ray reflRay = new Ray(hitPoint.x, hitPoint.y, hitPoint.z, reflection.x, reflection.y, reflection.z);
    Hit reflHit =  raySceneInteraction(reflRay, obj);


    PVector reflDir = new PVector();


    for (int j = 0; j<lights.size(); j++) {

      Light lightSource = lights.get(j);


      lightPos = lightSource.lightPosition;

      PVector.sub(lightPos, hitPoint, lightVec);
      lightVec.normalize();


      Ray shadRay = new Ray(hitPoint.x, hitPoint.y, hitPoint.z, lightVec.x, lightVec.y, lightVec.z);
      Hit shadHit = raySceneInteraction(shadRay, obj);

      if (isVisible(hit, lightSource)) {

        visibleLi = 1;

        PVector shadColor  = new PVector();
        shadColor.x = objects.get(obj).diffuseColor.x * lightSource.lightColor.x * max(0, lightVec.dot(shadHit.normal));
        shadColor.y = objects.get(obj).diffuseColor.y * lightSource.lightColor.y * max(0, lightVec.dot(shadHit.normal));
        shadColor.z = objects.get(obj).diffuseColor.z * lightSource.lightColor.z * max(0, lightVec.dot(shadHit.normal));

        float nDotL = hit.normal.dot(lightVec);
        PVector norm2 = PVector.mult(hit.normal, 2);
        PVector normNdotL = PVector.mult(norm2, nDotL);

        PVector.sub(normNdotL, lightVec, reflDir);
        reflDir.normalize();

        phongHL = pow(reflDir.dot(eyeRay), objects.get(obj).P);




        lightColorX += visibleLi*lightSource.lightColor.x*max(0, hit.normal.dot(lightVec))*objects.get(obj).diffuseColor.x + 
          lightSource.lightColor.x*objects.get(obj).specularColor.x*max(0, phongHL) + objects.get(obj).Krefl/**shadeHit(reflHit,obj, depth).x*/;


        lightColorY += visibleLi*lightSource.lightColor.y*max(0, hit.normal.dot(lightVec))*objects.get(obj).diffuseColor.y + 
          lightSource.lightColor.y*objects.get(obj).specularColor.y*max(0, phongHL) + objects.get(obj).Krefl/**shadeHit(reflHit,obj, depth).y*/;


        lightColorZ += visibleLi*lightSource.lightColor.z*max(0, hit.normal.dot(lightVec))*objects.get(obj).diffuseColor.z + 
          lightSource.lightColor.z*objects.get(obj).specularColor.z*max(0, phongHL) + objects.get(obj).Krefl/**shadeHit(reflHit,obj,depth).z*/;


        hit.hitColor.x = lightColorX;
        hit.hitColor.y = lightColorY;
        hit.hitColor.z = lightColorZ;
      }
    }
  }
  hit.hitColor.x += objects.get(obj).Car;
  hit.hitColor.y += objects.get(obj).Cag;
  hit.hitColor.z += objects.get(obj).Cab;

  depth++;

  return hit.hitColor;
}

boolean isVisible(Hit h, Light l) {
  boolean isVis = true;

  Ray visibleCheck = new Ray(h.hitPoint.x, h.hitPoint.y, h.hitPoint.z, l.lightPosition.x-h.hitPoint.x, l.lightPosition.y-h.hitPoint.y, l.lightPosition.z-h.hitPoint.z);
  if (visibleCheck.direction.dot(h.normal) < 0) {
    isVis = false;
  }
  return isVis;
}