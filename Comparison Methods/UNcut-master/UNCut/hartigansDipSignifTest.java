package graphClustering;

public class hartigansDipSignifTest {
	int[] triple;
	double dip,pvalue;

	public hartigansDipSignifTest(double[] x) throws Exception{
		int debug=0, min_is_0=1, n=x.length,nboot=1000,i,j;
		int[] triple = new int[3];
		int[] tmptriple = new int[3];
		double dip,pvalue;

		dip=dipStatistics.diptest(x,triple,debug,min_is_0);

		double[] bootDip = new double[nboot];
		double[] unifpdfboot = new double[n];


		for (i=0;i<nboot;i++){
			//Random generator = new Random();

			for (j=0;j<n;j++){
				//unifpdfboot[j]=generator.nextDouble();
				unifpdfboot[j]=Math.random();
			}
			myQuickSort.quickSort(unifpdfboot, 0, n-1);

			bootDip[i]=dipStatistics.diptest(unifpdfboot,tmptriple,debug,min_is_0);
		}

		myQuickSort.quickSort(bootDip, 0, nboot-1);
		double cnt=0;
		for (i=0;i<nboot;i++){
			if (dip<bootDip[i])
				cnt=cnt+1.0;
		}
		pvalue=cnt/nboot;



		this.triple=triple;
		this.dip=dip;
		this.pvalue=pvalue;
	}

	public int[] getTriple() {
		return triple;
	}

	public double getDip() {
		return dip;
	}

	public double getpValue() {
		return pvalue;
	}

}
