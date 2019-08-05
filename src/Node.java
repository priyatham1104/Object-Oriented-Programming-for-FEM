package FEM;

import iceb.jnumerics.MatrixFormat;
import iceb.jnumerics.Vector3D;


public class Node 
{
	private int[] dofNumbers = new int[3];			// is the number corresponding to the rows in the reduced stiffness matrix
	private Constraint c;	// for the nodes that doesn't have the constraints defined
	private Vector3D node_position;
	private Vector3D node_displacement;
	private Force f;
		
	public Node (double x1, double x2, double x3)
	{
		double[] position = new double [3];
		position[0] = x1;
		position[1] = x2;
		position[2] = x3;
		this.node_position = new Vector3D(position);
		Force f1 = new Force(0,0,0);
		this.setForce(f1);
		Constraint cons = new Constraint(true, true, true);
		this.setConstraint(cons);
		this.setDisplacement(0,0,0);
	}
	
	public void setConstraint(Constraint C)
	{
		this.c = C;
	}
	
	public Constraint getConstraint()
	{
		return this.c;
	}
	
	public void setForce(Force F)
	{
		this.f = F;
	}
	
	public Force getForce()
	{
		return this.f;
	}
	
	public int enumarateDOFs (int start)
	{
		for(int i =0 ; i < 3 ; i++)
		{
			if(this.c.isFree(i) == false)
			{
				dofNumbers[i] = -1;
			}
			else
			{
				dofNumbers[i] = start;
				start++;
			}
		}
		return start;
	}
	
	public int[] getDOFNumbers()
	{
		return dofNumbers;
	}
	
	public Vector3D getPosition()
	{
		return this.node_position;
	}
	
	public void setDisplacement(double u1, double u2, double u3)
	{
		double[] disp = new double [3];
		disp[0] = u1;
		disp[1] = u2;
		disp[2] = u3;
		this.node_displacement = new Vector3D(disp);
	}
	
	public void setposition(double u1, double u2, double u3)
	{
		double[] disp = new double [3];
		disp[0] = u1;
		disp[1] = u2;
		disp[2] = u3;
		this.node_position = new Vector3D(disp);
	}
	
	
	public Vector3D getDispacement()
	{
		return this.node_displacement;
	}
	
	public void print()
	{
		System.out.print(MatrixFormat.format(this.node_position));
	}
}
