class Light{
  PVector lightColor, lightPosition;
  float lightR, lightG, lightB, lightX, lightY, lightZ;
  
  Light(float lightR, float lightG, float lightB, float lightX, float lightY, float lightZ){
  
    this.lightR = lightR;
    this.lightG = lightG;
    this.lightB = lightB;
    
    this.lightX = lightX;
    this.lightY = lightY;
    this.lightZ = lightZ;
    
    lightColor = new PVector(lightR, lightG, lightB);
    lightPosition = new PVector(lightX, lightY, lightZ);
  }

}