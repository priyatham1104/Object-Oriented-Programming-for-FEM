package FEM;

import inf.v3d.obj.CylinderSet;
import inf.v3d.view.Viewer;
import iceb.jnumerics.*;
import iceb.jnumerics.HBReader.Matrix;
import iceb.jnumerics.lse.*;
import inf.text.ArrayFormat;
import inf.v3d.obj.Arrow;
import inf.v3d.obj.PolygonSet;

public class Visulaizer {

    private Structure structure;
    private Viewer viewer;
    private double scalingfactor_linear;
    public double scalingfactor_non_linear;
    private CylinderSet C;
  
    

    public Visulaizer(Structure S){
        this.structure = S;
        this.viewer = new Viewer();
    }

    public void drawElements(){
        CylinderSet C = new CylinderSet();
        for(int i = 0; i<this.structure.getNumberOfElements();i++){
            double [] p1 = this.structure.getElement(i).getNode1().getPosition().toArray();
            double [] p2 = this.structure.getElement(i).getNode2().getPosition().toArray();
            C.addCylinder(p1, p2, Math.sqrt((this.structure.getElement(i).getArea()/Math.PI)));
    }
        C.setRandomColor();
        this.viewer.addObject3D(C);
        this.viewer.setVisible(true);
    }
    
    public void drawelementforces(){
    	
    	double scalingfactor =0.0000001;
    	double [] k = new double[3];
    	k[0] = 0;
    	k[1] = 0;
    	k[2] = 1;
    	Vector3D n = new Vector3D(k);
    	{
    	for(int i = 0; i<this.structure.getNumberOfElements();i++){
    		
    		PolygonSet P = new PolygonSet();
    		Vector3D p1 = this.structure.getElement(i).getNode1().getPosition();
    		Vector3D p2 = this.structure.getElement(i).getNode2().getPosition();
    		double norm = (this.structure.getElement(i).getNode2().getPosition().subtract(this.structure.getElement(i).getNode1().getPosition())).normTwo();
            Vector3D d = (this.structure.getElement(i).getNode2().getPosition().subtract(this.structure.getElement(i).getNode1().getPosition())).multiply(1/norm);
            Vector3D p = n.vectorProduct(d);
            Vector3D s1 = p1.add((p.multiply(scalingfactor*(this.structure.getElement(i).computeForce()))));
            Vector3D s2 = p2.add((p.multiply(scalingfactor*(this.structure.getElement(i).computeForce()))));
          P.insertVertex(p1.toArray()[0],p1.toArray()[1],p1.toArray()[2],1);
          P.insertVertex(p2.toArray()[0],p2.toArray()[1],p2.toArray()[2],1);
          P.insertVertex(s2.toArray()[0],s2.toArray()[1],s2.toArray()[2],1);
          P.insertVertex(s1.toArray()[0],s1.toArray()[1],s1.toArray()[2],1); 
            P.polygonComplete();
            P.setRandomColor();
            this.viewer.addObject3D(P);
    }
    	
    	
         this.viewer.setVisible(true);
    	}
    	
    }

    public void Elementforces(){
        for(int i = 0; i<this.structure.getNumberOfElements();i++){
            double [] p1 = this.structure.getElement(i).getNode1().getPosition().toArray();
            double [] p2 = this.structure.getElement(i).getNode2().getPosition().toArray();
            Arrow p = new Arrow(p1,p2);
            p.setRandomColor();
            p.setRadius(0.01);
            this.viewer.addObject3D(p);
        }
        this.viewer.setVisible(true);
    }
    
    public void setscalingfactor_linear(double k){
    	 this.scalingfactor_linear = k;
    	
    	
    }
    
    public void setscalingfactor_non_linear(double k){
   	 this.scalingfactor_non_linear = k;
   	
   	
   }
    public void Constraints(){
    	for(int i=0; i<this.structure.getNumberOfNodes(); i++){
    		for(int j=0; j<3; j++){
    			if (!this.structure.getNode(i).getConstraint().isFree(j)){
    				if(j==0){
    		            double [] p1 = this.structure.getNode(i).getPosition().toArray();
    		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
    		            p1[0] = p1[0]-0.5;
    		            Arrow p = new Arrow(p1,p2);
    		            p.setRadius(0.02);
    		            p.setColor("blue");
    		            this.viewer.addObject3D(p);
    				}
    				if(j==1){
    		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
    		            double [] p1 = this.structure.getNode(i).getPosition().toArray();
    		            p1[1] = p1[1]-0.5;
    		            Arrow p = new Arrow(p1,p2);
    		            p.setRadius(0.02);
    		            p.setColor("green");
    		            this.viewer.addObject3D(p);
    				}
    				if(j==2){
    		            double [] p1 = this.structure.getNode(i).getPosition().toArray();
    		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
    		            p1[2] = p1[2]-0.5;
    		            Arrow p = new Arrow(p1,p2);
    		            p.setRadius(0.02);
    		            p.setColor("red");
    		            this.viewer.addObject3D(p);
    				}
    			}
    		}
    		this.viewer.setVisible(true);
    	}
    }
    
    public void drawNodalForce(){
    	for(int i=0; i<this.structure.getNumberOfNodes(); i++){
    		for(int j=0; j<3; j++){
    			if(Math.abs((this.structure.getNode(i).getForce().getComponent(j)))>0){
    				
    				if(j==0){
    					if(this.structure.getNode(i).getForce().getComponent(j)<0){
    						double [] p1 = this.structure.getNode(i).getPosition().toArray();
        		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
        		            p1[0] = p1[0]-0.5;
        		            Arrow p = new Arrow(p2,p1);
        		            p.setRadius(0.05);
        		            p.setColor("blue");
        		            this.viewer.addObject3D(p);
        				}
    					if(this.structure.getNode(i).getForce().getComponent(j)>0){
    						double [] p1 = this.structure.getNode(i).getPosition().toArray();
        		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
        		            p2[0] = p2[0]+0.5;
        		            Arrow p = new Arrow(p1,p2);
        		            p.setRadius(0.05);
        		            p.setColor("blue");
        		            this.viewer.addObject3D(p);
    						
    					}
    					}
    		            
    				if(j==1){
    					if(this.structure.getNode(i).getForce().getComponent(j)<0){
    						double [] p1 = this.structure.getNode(i).getPosition().toArray();
        		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
        		            p2[1] = p2[1]+0.5;
        		            Arrow p = new Arrow(p2,p1);
        		            p.setRadius(0.02);
        		            p.setColor("blue");
        		            this.viewer.addObject3D(p);
        				}
    					if(this.structure.getNode(i).getForce().getComponent(j)>0){
    						double [] p1 = this.structure.getNode(i).getPosition().toArray();
        		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
        		            p2[1] = p2[1]+0.5;
        		            Arrow p = new Arrow(p1,p2);
        		            p.setRadius(0.05);
        		            p.setColor("blue");
        		            this.viewer.addObject3D(p);
    						
    					}
    					}
    		            
    				if(j==2){
    					if(this.structure.getNode(i).getForce().getComponent(j)<0){
    						double [] p1 = this.structure.getNode(i).getPosition().toArray();
        		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
        		            p1[2] = p1[2]+0.5;
        		            Arrow p = new Arrow(p1,p2);
        		            p.setRadius(0.05);
        		            p.setColor("blue");
        		            this.viewer.addObject3D(p);
        				}
    					if(this.structure.getNode(i).getForce().getComponent(j)>0){
    						double [] p1 = this.structure.getNode(i).getPosition().toArray();
        		            double [] p2 = this.structure.getNode(i).getPosition().toArray();
        		            p1[2] = p1[2]+0.5;
        		            Arrow p = new Arrow(p1,p2);
        		            p.setRadius(0.05);
        		            p.setColor("blue");
        		            this.viewer.addObject3D(p);
    						
    					}
    					}
    			}
    		}
    	}
    	this.viewer.setVisible(true);
    	
    }
    
    public void setscalingfactorauto(){
    	System.out.print(ArrayFormat.format(this.structure.getU()));
    	double scale = 0;
    	for(int i=0; i<this.structure.getU().length; i++){
    		scale = Math.max(scale, Math.abs(this.structure.getU()[i]));
    	}
    	System.out.print("Scaling factor is"+ scale);
    	this.scalingfactor_linear = 1/(scale*5);
    	this.scalingfactor_non_linear = 1/(scale*5);
    	
    }
    
    public void drawdeformedconfiduration_linear(){
       
        for(int i=0;i<this.structure.getNumberOfNodes();i++){
            double [] d = this.structure.getNode(i).getDispacement().toArray();
            double [] p = this.structure.getNode(i).getPosition().toArray();
            for(int j=0;j<3;j++){
                p[j] = p[j] + this.scalingfactor_linear*d[j];
            }
            this.structure.getNode(i).setposition(p[0], p[1], p[2]);

        }
        this.drawElements();
//        this.C.setColor("blue");
    }
    public void drawdeformedconfiduration_non_linear(){
        
        for(int i=0;i<this.structure.getNumberOfNodes();i++){
            double [] d = this.structure.getNode(i).getDispacement().toArray();
            double [] p = this.structure.getNode(i).getPosition().toArray();
            for(int j=0;j<3;j++){
                p[j] = p[j] + this.scalingfactor_non_linear*d[j];
            }
            this.structure.getNode(i).setposition(p[0], p[1], p[2]);

        }
        this.drawElements();
//        this.C.setColor("blue");
    }
   
}