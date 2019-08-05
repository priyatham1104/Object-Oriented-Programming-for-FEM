package FEM;

import iceb.jnumerics.*;
//import iceb.jnumerics.IMatrix;
//import inf.jlinalg.IMatrix;
import inf.text.ArrayFormat;

public class Element 
{
	
	private double a;
	private double eM;
	private int[] dofNumbers = new int[6];
	private Node n1 , n2;
	
	public Element (double e, double a, Node n1, Node n2)
	{
		this.eM = e;
		this.a = a;
		this.n1 = n1;
		this.n2 = n2;
	}
	
	public IMatrix computeStiffnessMatrix(){
	    Vector3D E1 = this.n1.getPosition().subtract(this.n2.getPosition());
	    Vector3D E2 = this.n2.getPosition().subtract(this.n1.getPosition());
	    IMatrix a = E1.dyadicProduct(E1);
	    IMatrix b = E1.dyadicProduct(E2);
	    IMatrix c = E2.dyadicProduct(E1);
	    IMatrix d = E2.dyadicProduct(E2);
	    IMatrix K = new Array2DMatrix(6,6);
	    double s = this.eM*this.a/(Math.pow(this.getLength(),3));
	for(int i=0; i < 3 ;i++){
	    for(int j=0;j<3;j++){
	        K.set(i, j, s*a.get(i, j));
	    }
	}
	for(int i=0; i < 3 ;i++){
	    for(int j=3;j<6;j++){
	        K.set(i, j, s*b.get(i, j-3));
	    }
	}
	for(int i=3; i < 6 ;i++){
	        for(int j=0;j<3;j++){
	            K.set(i, j, s*c.get(i-3, j));
	        }
	}
	for(int i=3; i < 6 ;i++){
	    for(int j=3;j<6;j++){
	        K.set(i, j, s*d.get(i-3, j-3));
	    }
	}
	System.out.print(K);
	return K;
	}
	
	public double[] compute_R_internal(){
		Vector3D X1 = this.n1.getPosition().add(this.n1.getDispacement()) ;
		Vector3D X2 = this.n2.getPosition().add(this.n2.getDispacement()) ;
		 Vector3D E1 = X1.subtract(X2);
		    Vector3D E2 = X2.subtract(X1);
		double l = (X1.subtract(X2)).normTwo();
		double L = this.getLength();
		double k = this.getEModulus()*this.getArea()/(2*Math.pow(L,3))*(l*l - L*L);
        double [] R_internal = new double[6];
	
		R_internal[0] = k*E1.get(0);
		R_internal[1]= k*E1.get(1);
		R_internal[2]=k*E1.get(2);
		R_internal[3]= k*E2.get(0);
		R_internal[4]= k*E2.get(1);
		R_internal[5]= k*E2.get(2);
//		System.out.print("\nthe R_internal vector for element is" + "\n");
//		System.out.print(ArrayFormat.format(R_internal));
		return R_internal;
	}

	public IMatrix Compute_non_linear_StiffnessMatrix(){
		Vector3D X1 = this.n1.getPosition().add(this.n1.getDispacement()) ;
		Vector3D X2 = this.n2.getPosition().add(this.n2.getDispacement()) ;
	    Vector3D E1 = X1.subtract(X2);
	    Vector3D E2 = X2.subtract(X1);
	    IMatrix a = E1.dyadicProduct(E1);
	    IMatrix b = E1.dyadicProduct(E2);
	    IMatrix c = E2.dyadicProduct(E1);
	    IMatrix d = E2.dyadicProduct(E2);
	    IMatrix K_m = new Array2DMatrix(6,6);
	    IMatrix K = new Array2DMatrix(6,6);
	    double s = this.eM*this.a/(Math.pow(this.getLength(),3));
	for(int i=0; i < 3 ;i++){
	    for(int j=0;j<3;j++){
	        K_m.set(i, j, s*a.get(i, j));
	    }
	}
	for(int i=0; i < 3 ;i++){
	    for(int j=3;j<6;j++){
	        K_m.set(i, j, s*b.get(i, j-3));
	    }
	}
	for(int i=3; i < 6 ;i++){
	        for(int j=0;j<3;j++){
	            K_m.set(i, j, s*c.get(i-3, j));
	        }
	}
	for(int i=3; i < 6 ;i++){
	    for(int j=3;j<6;j++){
	        K_m.set(i, j, s*d.get(i-3, j-3));
	    }
	}
	
	IMatrix K_g = new Array2DMatrix(6,6);
	IMatrix A = new Array2DMatrix(6,6);
    IMatrix I = new Array2DMatrix(3,3);
    
    		for(int i=0;i<3;i++){
    			for(int j=0;j<3;j++){
    				if(i==j){
    	    			I.set(i, j, 1);	
    			}
    				else I.set(i, j, 0);
    		
    		}
    		}
    		
    		for(int i=0;i<6;i++){
    			for(int j=0;j<6;j++){
    				if(i<3 && j<3){
    					A.set(i, j, I.get(i, j));
    				}
    				if(i<3 && j>=3){
    					A.set(i, j, -I.get(i, j-3));
    				}
    				if(i>=3 && j<3){
    					A.set(i, j, -I.get(i-3, j));
    				}
    				if(i>=3 && j>=3){
    					A.set(i, j, I.get(i-3, j-3));
    				}
    			}
    		}
double l = (X1.subtract(X2)).normTwo();
double L = this.getLength();
double k = this.getEModulus()*this.getArea()/(2*Math.pow(L,3))*(l*l - L*L);

for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
		K_g.set(i, j, k*A.get(i, j));
	}	
	}

for(int i=0;i<6;i++){
	for(int j=0;j<6;j++){
		K.set(i, j, (K_m.get(i, j)+K_g.get(i, j)));
	}	
	}
//System.out.println("The Geometric Stiffness Matrix: "+"\n"+ K_g);
//System.out.println("The Matrial Stiffness Matrix: "+"\n"+ K_m);
//System.out.println("The Tangential Stiffness Matrix: "+"\n"+ K);
return K;

	}
	
	public IMatrix compute_non_linear_StiffnessMatrix(){
		Vector3D E1 = this.n1.getPosition().subtract(this.n2.getPosition());
	    Vector3D E2 = this.n2.getPosition().subtract(this.n1.getPosition());
	    IMatrix A = new Array2DMatrix(6,6);
	    IMatrix I = new Array2DMatrix(3,3);
	    
	    		for(int i=0;i<3;i++){
	    			for(int j=0;j<3;j++){
	    				if(i==j){
	    	    			I.set(i, j, 1);	
	    			}
	    				else I.set(i, j, 0);
	    		
	    		}
	    		}
	    		
	    		for(int i=0;i<6;i++){
	    			for(int j=0;j<6;j++){
	    				if(i<3 && j<3){
	    					A.set(i, j, I.get(i, j));
	    				}
	    				if(i<3 && j>=3){
	    					A.set(i, j, -I.get(i, j-3));
	    				}
	    				if(i>=3 && j<3){
	    					A.set(i, j, -I.get(i-3, j));
	    				}
	    				if(i>=3 && j>=3){
	    					A.set(i, j, I.get(i-3, j-3));
	    				}
	    			}
	    		}

	    		
	    		System.out.print(A);
	    		System.out.print(I);
	    	return A;	
	}
	
	
	
	public void enumerateDOFs()
	{
		for(int i = 0; i < 6; i++)
		{
			if(i<3)
			{
			dofNumbers[i] = n1.getDOFNumbers()[i];
			}
			else
			{
			dofNumbers[i] = n2.getDOFNumbers()[i-3];
			}
		}
	}
			
	
	public int[] getDOFNumbers()
	{
		return this.dofNumbers;
	}
	
	public Vector3D getE1()
	{
		Vector3D l;
		l = this.n2.getPosition().subtract(this.n1.getPosition());
		return l;
	}
	
	public double getLength()
	{
	
		double L;
		L = this.getE1().normTwo();
		return L;
	}
	
	public double computeForce()				
	{
		Vector3D disp_vect = this.getNode2().getPosition().subtract(this.getNode1().getPosition());
		double u1 = this.getNode1().getDispacement().scalarProduct(disp_vect)/this.getLength();
		double u2 = this.getNode2().getDispacement().scalarProduct(disp_vect)/this.getLength();
		double F = (this.eM*this.getArea()/this.getLength())*(u2-u1);
	
		return F;
	}
	
	public Node getNode1()
	{
		return this.n1;
	}
	
	public Node getNode2()
	{
		return this.n2;
	}
	
	public double getArea()
	{
		return this.a;
	}
	
	public double getEModulus()
	{
		return this.eM;
	}
	
	public void print()
	{
		System.out.print(this.getEModulus());
		System.out.print("\t");
		System.out.print(this.getArea());
		System.out.print("\t");
		System.out.print(this.getLength());
	}
}
