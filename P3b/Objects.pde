abstract class Objects {
  //PVector objectColor;
  PVector diffuseColor, specularColor;
  float objectID, Cdr, Cdg, Cdb, Csr, Csg, Csb, P, Krefl, Car, Cag, Cab;

  Objects(float Cdr, float Cdg, float Cdb, float Csr, float Csg, float Csb, float P, float Krefl, float Car, float Cag, float Cab) {
    objectID = 0;
    this.Cdr = Cdr;
    this.Cdg = Cdg;
    this.Cdb = Cdb;
    this.Csr = Csr;
    this.Csg = Csg;
    this.Csb = Csb;
    this.Krefl = Krefl;
    this.P = P;
    this.Car = Car;
    this.Cag = Cag;
    this.Cab = Cab;
 
    diffuseColor = new PVector(Cdr, Cdg, Cdb);
    specularColor = new PVector(Csr, Csg, Csb);

  }
}

class Sphere extends Objects {
  PVector sphereCenter;
  float radius;

  Sphere(float Cdr, float Cdg, float Cdb, float Csr, float Csg, float Csb, float P, float Krefl, float Car, float Cag, float Cab) {
    super(Cdr, Cdg, Cdb, Csr, Csg, Csb, P, Krefl, Car, Cag, Cab);

    objectID = 1;

    sphereCenter = new PVector(sphereX, sphereY, sphereZ);

    radius = sphereR;
  }
}

class Triangle extends Objects {
  PVector p1, p2, p3;
  
  Triangle(float Cdr, float Cdg, float Cdb, float Csr, float Csg, float Csb, float P, float Krefl, float Car, float Cag, float Cab) {
    super(Cdr, Cdg, Cdb, Csr, Csg, Csb, P, Krefl, Car, Cag, Cab);
    objectID = 2;

    p1 = new PVector(0,0,0);
    p2 = new PVector(0, 0, 0);
    p3 = new PVector(0,0,0);
    
    
  }
  
}