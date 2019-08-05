package FEM;

import iceb.jnumerics.Array2DMatrix;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.QuadraticMatrixInfo;
import iceb.jnumerics.SolveFailedException;
import iceb.jnumerics.lse.GeneralMatrixLSESolver;
import iceb.jnumerics.lse.ILSESolver;
import inf.text.ArrayFormat;

import java.util.ArrayList;

public class Structure {

	private Node n;
	private Element E;
	private ArrayList<Node> nodes = new ArrayList<Node>();
	private ArrayList<Element> elements = new ArrayList<Element>();
	private double[] U;

	public Structure() {
	}
	
	public double[] getU(){
		return this.U;
	}

	public Node addNode(double x1, double x2, double x3) {
		nodes.add(this.n = new Node(x1, x2, x3));
		return this.n;
	}

	public Element addElement(double e, double a, int n1, int n2) {
		elements.add(this.E = new Element(e, a, nodes.get(n1), nodes.get(n2)));
		return this.E;
	}

	public int getNumberOfNodes() {
		return nodes.size();
	}

	public Node getNode(int i) {
		return this.nodes.get(i);
	}

	public int getNumberOfElements() {
		return this.elements.size();
	}

	public Element getElement(int i) {
		return this.elements.get(i);
	}

//	public void printStructure() {
//		System.out.print("The number of elements in the structure = "
//				+ elements.size());
//		System.out.print("\nThe Elements of the structure are: \n");
//		for (int i = 0; i < this.elements.size(); i++) {
//			this.elements.get(i).print();
//			System.out.print("\n");
//		}
//		System.out.print("The number of nodes in the structure = "
//				+ nodes.size());
//		System.out.print("\nThe Nodes of the structure are: \n");
//		for (int i = 0; i < this.nodes.size(); i++) {
//			this.nodes.get(i).print();
//			System.out.print("\n");
//		}
//	}

	private int enumerateDOFs() {
		int i = 0;
//		System.out.print("\nThe Node degrees of freedom: \n");
		for (int j = 0; j < this.nodes.size(); j++) {
			int p = this.nodes.get(j).enumarateDOFs(i);
			int[] dof = this.nodes.get(j).getDOFNumbers();
//			System.out.print(ArrayFormat.format(dof));
//			System.out.print("\n");
			i = p;
		}
		int n = 0;
//		System.out.print("\nThe Element degrees of freedom: \n");
		for (int k = 0; k < this.elements.size(); k++) {
			this.elements.get(k).enumerateDOFs();
			int[] dof = this.elements.get(k).getDOFNumbers();
//			System.out.print(ArrayFormat.format(dof));
//			System.out.print("\n");
		}
		return i;
	}
	
	public void Listing() {
		System.out.print("\n" + "\n" + "Listing Structure" + "\n");
		System.out.print("\n"+"The number of nodes in the structure = "
				+ this.getNumberOfNodes() + "\n");
		System.out.print("\n"+"Nodes \n");
		System.out.print("The number of elements in the structure = "
				+ this.getNumberOfElements() + "\n");
		System.out.print("Nodes \n");
		System.out.print("idx	x1		x2 		x3 \n");
		for (int i = 0; i < this.nodes.size(); i++) {
			System.out.print(i);
			System.out.print(this.nodes.get(i).getPosition() + "\n");
		}
		
		System.out.print("\n"+"Constrains \n");
		System.out.print("node	u1	u2 	u3 ");
		for (int i = 0; i < this.nodes.size(); i++) {
			System.out.print("\n" + i + "\t");
			this.nodes.get(i).getConstraint().print();
			
		}
		System.out.print("\n");
		
		System.out.print("\n"+"Forces \n");
		System.out.print("node	r1	r2 	r3 \n");
		
		for (int i = 0; i < this.nodes.size(); i++) {
			if(!(this.nodes.get(i).getForce().getComponent(0) == 0 &&this.nodes.get(i).getForce().getComponent(1) == 0&&this.nodes.get(i).getForce().getComponent(2) == 0 )){
			System.out.print(i + "\t");
			System.out.print(this.nodes.get(i).getForce().getComponent(0) + "\t");
			System.out.print(this.nodes.get(i).getForce().getComponent(1) + "\t");
			System.out.print(this.nodes.get(i).getForce().getComponent(2) + "\t" + "\n");
			}
		}
		
		
		System.out.print("\n"+"Elements \n");
		System.out.print("idx	E	A 		Length \n");
		for (int i = 0; i < this.elements.size(); i++) {
			System.out.print(i);
			System.out.print("\t" + this.elements.get(i).getEModulus()	+ "\t" + this.elements.get(i).getArea() + "\t" + this.elements.get(i).getLength() + "\n");
		}
		
		System.out.print("\n"+"Listing Analysis Results \n");
		System.out.print("\n"+"Displacements \n");
		System.out.print("node		u1		u2 		u3 \n");
		for (int i = 0; i < this.nodes.size(); i++) {
			System.out.print(i);
			System.out.print("\t" + this.nodes.get(i).getDispacement()	+ "\n");
		}
		
		System.out.print("Element Forces \n");
		System.out.print("elem		forces \n");
		for (int i = 0; i < this.elements.size(); i++) {
			System.out.print(i);
			System.out.print("\t" + this.elements.get(i).computeForce()	+ "\n");
		}
	}
	
	

	private void assembleLoadVector(double[] rGlobal) {
		int k = 0;
		for (int i = 0; i < this.nodes.size(); i++) {
			for (int j = 0; j < 3; j++) {
				if (this.nodes.get(i).getConstraint().isFree(j)) {
					rGlobal[k] = this.nodes.get(i).getForce().getComponent(j);
					k++;
				}
			}
		}
//		System.out.print("\nThe assembled Load Vector is: \n");
//		System.out.print(ArrayFormat.format(rGlobal));
//		System.out.print("\n");
	}

	private void assembleStiffnessMatrix(IMatrix kGlobal) {
		for (int i = 0; i < this.elements.size(); i++) {
			this.elements.get(i).enumerateDOFs();
			for (int r = 0; r < 6; r++) {
				for (int c = 0; c < 6; c++) {
					if (this.elements.get(i).getDOFNumbers()[r] >= 0) {
						if (this.elements.get(i).getDOFNumbers()[c] >= 0) {
							kGlobal.set(
									this.elements.get(i).getDOFNumbers()[r],
									this.elements.get(i).getDOFNumbers()[c],
									kGlobal.get(this.elements.get(i)
											.getDOFNumbers()[r], this.elements
											.get(i).getDOFNumbers()[c])
											+ this.elements.get(i)
													.computeStiffnessMatrix()
													.get(r, c));
						}
					}
				}
			}
		}
		System.out.print("\nThe assembled Stiffness Matrix is: \n");
		System.out.print(kGlobal);
		System.out.print("\n");
	}

	private double [] assembleInternalLoadvector() {
		double [] rinternal = new double [this.enumerateDOFs()];
//		for (int k=0;k<this.enumerateDOFs();k++){
//			rinternal[k]=0;
//		}
		
		for (int i = 0; i < this.elements.size(); i++) {
			this.elements.get(i).enumerateDOFs();
			for (int r = 0; r < 6; r++) {
					if (this.elements.get(i).getDOFNumbers()[r] >= 0) {						
									
									rinternal[this.elements.get(i).getDOFNumbers()[r]] = rinternal[this.elements.get(i).getDOFNumbers()[r]] + this.elements.get(i).compute_R_internal()[r];
						}
					}
				
		}
//		System.out.print("\nThe internal load vector is: \n");
//		System.out.print(ArrayFormat.format(rinternal));
//		System.out.print("\n");
		return rinternal;
	}

	private double[] luSolver(IMatrix kGlobal, double[] rGlobal) {
		int m = rGlobal.length;
		ILSESolver solver = new GeneralMatrixLSESolver();
		QuadraticMatrixInfo kInfo = solver.getAInfo();

		IMatrix K = solver.getA();
		kInfo.setSize(m);
		solver.initialize();

		for (int i = 0; i < m; i++) {

			for (int j = 0; j < m; j++) {
				K.set(i, j, kGlobal.get(i, j));
			}
		}
		double[] R = rGlobal;

		try {
			solver.solve(R);
		} catch (SolveFailedException e) {
			System.out.println("Solver failed: " + e.getMessage());
		}

//		System.out.println("Solution X");
//		System.out.print(ArrayFormat.format(R));
		return R;
	}

	public void solve()
    {
            int m = this.enumerateDOFs();
            double[] rGlobal = new double[m];
            this.assembleLoadVector(rGlobal);
            IMatrix kGlobal = new Array2DMatrix(m,m);
            this.assembleStiffnessMatrix(kGlobal);
            this.U = this.luSolver(kGlobal, rGlobal);
             for(int i =0,k=0;i<this.nodes.size();i++)
             {
                    double [] q = new double[3];
                    for(int j=0;j<3;j++){
                    if(this.nodes.get(i).getDOFNumbers()[j] >= 0){
                        q[j] = this.U[k];
                        k = k+1;
                    }
                    this.nodes.get(i).setDisplacement(q[0],q[1], q[2]);
                    }
               }

            System.out.print("\nThe assembled Load Vector is: \n");
            System.out.print(ArrayFormat.format(rGlobal));
            System.out.print("\n\nThe assembled Stiffness Matrix is: \n");
            System.out.print(kGlobal);
	}
	
public double Arraynorm(double [] k){
	double p = 0;
	for(int i=0;i<k.length;i++){
		p = p + k[i]*k[i];

	}
	return Math.sqrt(p);
}
	
	private void Assemble_Non_Linear_StiffnessMatrix(IMatrix kGlobal_non_linear) {
		for (int i = 0; i < this.elements.size(); i++) {
			this.elements.get(i).enumerateDOFs();
			for (int r = 0; r < 6; r++) {
				for (int c = 0; c < 6; c++) {
					if (this.elements.get(i).getDOFNumbers()[r] >= 0) {
						if (this.elements.get(i).getDOFNumbers()[c] >= 0) {
							kGlobal_non_linear.set(
									this.elements.get(i).getDOFNumbers()[r],
									this.elements.get(i).getDOFNumbers()[c],
									kGlobal_non_linear.get(this.elements.get(i)
											.getDOFNumbers()[r], this.elements
											.get(i).getDOFNumbers()[c])
											+ this.elements.get(i)
													.Compute_non_linear_StiffnessMatrix()
													.get(r, c));
						}
					}
				}
			}
		}
//		System.out.print("\nThe assembled Non Linear Stiffness Matrix is: \n");
//		System.out.print(kGlobal_non_linear);
//		System.out.print("\n");
	}



private double[] luSolver_non_linear(IMatrix kGlobal_non_linear, double[] rGlobal) {
	int m = rGlobal.length;
	ILSESolver solver = new GeneralMatrixLSESolver();
	QuadraticMatrixInfo kInfo = solver.getAInfo();

	IMatrix K = solver.getA();
	kInfo.setSize(m);
	solver.initialize();

	for (int i = 0; i < m; i++) {

		for (int j = 0; j < m; j++) {
			K.set(i, j, kGlobal_non_linear.get(i, j));
		}
	}
	double[] R = rGlobal;

	try {
		solver.solve(R);
	} catch (SolveFailedException e) {
		System.out.println("Solver failed: " + e.getMessage());
	}

	System.out.println("Solution X");
	System.out.print(ArrayFormat.format(R));
	return R;
}


public void solve_non_linear()
{
	double b =0;
	double load_increment = 0.2;
	this.U = new double[this.enumerateDOFs()];
	double [] U_lambda = new double[this.enumerateDOFs()];
	double [] U_old = new double[this.enumerateDOFs()];
    double [] U_new = new double[this.enumerateDOFs()];
    System.out.print("\nThe initial displacement vector for  iteration"+"\n");
    System.out.print(ArrayFormat.format(U_old));
	double lambda =0;
	
	lambda = lambda + load_increment;
	for(;lambda<=1;){
		System.out.print("\nThe lambda for the iteration is"+"\n");		
System.out.print(lambda);
double ETA =100;
	while(ETA>1E-05){
		b = b+1;
        int m = this.enumerateDOFs();     
        double[] r_ext = new double[m];
        double[] r_int = new double[m];
        double[] r_global = new double[m];
        r_int = this.assembleInternalLoadvector();
      System.out.print("\nThe assembled internal Load Vector for iteration"+" "+b+"\n");
      System.out.print(ArrayFormat.format(r_int));
        this.assembleLoadVector(r_ext);
        System.out.print("\nThe assembled external Load Vector for iteration"+" "+b+"\n");
        System.out.print(ArrayFormat.format(r_ext));
        for(int a=0;a<r_ext.length;a++){
        	r_global[a] = lambda*r_ext[a]-this.assembleInternalLoadvector()[a];
        }
        System.out.print("\nThe assembled global Load Vector for iteration"+" "+b+"\n");
        System.out.print(ArrayFormat.format(r_global));
        IMatrix kGlobal_non_linear = new Array2DMatrix(m,m);
        this.Assemble_Non_Linear_StiffnessMatrix(kGlobal_non_linear);
//        System.out.print("\nThe assembled Stiffness Matrix for iteration"+" "+b+"\n");
//        System.out.print(kGlobal_non_linear); 
       double [] delta_U = this.luSolver(kGlobal_non_linear,r_global);
       System.out.print("\nThe delta_u for iteration"+" "+b+"\n");
       System.out.print(ArrayFormat.format(delta_U)); 
        for(int n=0;n<m;n++){
        	U_new[n] = U_old[n]+delta_U[n];
        	
        } 
        double [] delta_U_lambda = new double[this.enumerateDOFs()];
        for(int g=0;g<this.enumerateDOFs();g++){
        	delta_U_lambda[g] = U_new[g]-U_lambda[g];
        }
        ETA = this.Arraynorm(delta_U)/this.Arraynorm(delta_U_lambda);
        System.out.print("\n ETA"+"\n");
        System.out.print(ETA);
        System.out.print("\nThe U_old for iteration"+" "+b+"\n");
        System.out.print(ArrayFormat.format(U_old));     
        System.out.print("\nThe U_new for iteration"+" "+b+"\n");
        System.out.print(ArrayFormat.format(U_new));
        for(int n=0;n<m;n++){
        	U_old[n] = U_new[n];	
        }

         for(int i =0,k=0;i<this.nodes.size();i++)
         {
                double [] q = new double[3];
                for(int j=0;j<3;j++){
                if(this.nodes.get(i).getDOFNumbers()[j] >= 0){
                    q[j] = U_new[k];
                    k = k+1;
                }
                this.nodes.get(i).setDisplacement(q[0],q[1], q[2]);
                }
                

                
//                for(int p=0, l =0;p<U_new.length;p++){
//
//                }
//                this.Assemble_Non_Linear_StiffnessMatrix(kGlobal_non_linear);
//            	System.out.print("\n\nThe assembled Stiffness Matrix for iteration no:" + n +"is\n");
//            	 System.out.print(kGlobal_non_linear);
           }
         System.out.print("\n"+"Listing Analysis Results \n");
 		System.out.print("\n"+"Displacements \n");
 		System.out.print("node		u1		u2 		u3 \n");
 		for (int l = 0; l < this.nodes.size(); l++) {
 			System.out.print(l);
 			System.out.print("\t" + this.nodes.get(l).getDispacement()	+ "\n");
 		}
        }
	for(int q=0;q<this.enumerateDOFs();q++){
		U_lambda[q] = U_new[q] ;
	}
	  System.out.print("\nThe U_lambda for iteration"+"\n");
      System.out.print(ArrayFormat.format(U_lambda));
	lambda = lambda + load_increment;
        }
	
	} } 



