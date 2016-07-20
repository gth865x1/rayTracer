class Ray{
  PVector origin, direction;
  Ray(){
    origin = new PVector(0, 0, 0);
    
    direction = new PVector(0, 0, -1);

  }
  
  Ray(float dx, float dy, float dz){
    origin = new PVector(0, 0, 0);

    direction = new PVector(dx-origin.x, dy-origin.y, dz-origin.z);
    direction.normalize();
 
  
  }
  
  Ray(float origX, float origY, float origZ, float dx, float dy, float dz){
    origin = new PVector(origX, origY, origZ);
    
    direction = new PVector(dx-origin.x, dy-origin.y, dz-origin.z);
    direction.normalize();
     
  }




}