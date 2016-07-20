class Hit{
  PVector hitPoint;
  PVector hitColor;
  int hitObject;
  float time;
  Ray eyeRay;
  PVector normal;
  
  Hit(){
    
  }
  
  Hit(float t, Ray r, int hitObNum, PVector hitNormal){
        time = t;
    
    eyeRay = r;
    
    normal = hitNormal;
    
    float x0 = r.origin.x;
    float y0 = r.origin.y;
    float z0 = r.origin.z;
    
    float tDx = time*r.direction.x;
    float tDy = time*r.direction.y;
    float tDz = time*r.direction.z;


    
    float hitX = x0+tDx;
    float hitY = y0+tDy;
    float hitZ = z0+tDz;
    hitPoint = new PVector(hitX, hitY, hitZ);
    
    hitColor = new PVector(0, 0, 0);
    
    hitObject = hitObNum;
  }

}