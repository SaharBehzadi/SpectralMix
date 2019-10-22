package graphClustering;

public class myQuickSort {
	//	public static void main(String[] args) {
	//
	//        //int[] a = { 1, 23, 45, 2, 8, 134, 9, 4, 2000 };
	//        int a[][]={{23,44,1,2009,2,88,123,7,999,1040,88},{1,2,3,4,5,6,7,8,9,10,11}};
	//        quickSort(a, 0, a[0].length - 1);
	//        for (int i=0;i<2;i++){
	//        	for (int j=0;j<11;j++){
	//        		System.out.print(a[i][j]);
	//        		System.out.print(' ');
	//        	}
	//        	System.out.print('\n');
	//        }
	//        
	//    }
	public static void quickSort(int[][] a, int p, int r)
	{
		if(p<r)
		{
			int q=partition(a,p,r);
			quickSort(a,p,q);
			quickSort(a,q+1,r);
		}
	}

	private static int partition(int[][] a, int p, int r) {

		int x = a[0][p];
		int i = p-1 ;
		int j = r+1 ;

		while (true) {
			i++;
			while ( i< r && a[0][i] < x)
				i++;
			j--;
			while (j>p && a[0][j] > x)
				j--;

			if (i < j)
				swap(a, i, j);
			else
				return j;
		}
	}

	private static void swap(int[][] a, int i, int j) {
		// TODO Auto-generated method stub
		int row=a.length;
		int temp;
		for (int i1=0;i1<row;i1++){
			temp = a[i1][i];
			a[i1][i] = a[i1][j];
			a[i1][j] = temp;
		}
	}
	
	//===========================================================
	public static void quickSort(double[][] a, int p, int r)
	{
		if(p<r)
		{
			int q=partition(a,p,r);
			quickSort(a,p,q);
			quickSort(a,q+1,r);
		}
	}

	private static int partition(double[][] a, int p, int r) {

		double x = a[0][p];
		int i = p-1 ;
		int j = r+1 ;

		while (true) {
			i++;
			while ( i< r && a[0][i] < x)
				i++;
			j--;
			while (j>p && a[0][j] > x)
				j--;

			if (i < j)
				swap(a, i, j);
			else
				return j;
		}
	}

	private static void swap(double[][] a, int i, int j) {
		// TODO Auto-generated method stub
		int row=a.length;
		double temp;
		for (int i1=0;i1<row;i1++){
			temp = a[i1][i];
			a[i1][i] = a[i1][j];
			a[i1][j] = temp;
		}
	}
	//=============================================================

	public static void quickSort(double[] a, int p, int r)
	{
		if(p<r)
		{
			int q=partition(a,p,r);
			quickSort(a,p,q);
			quickSort(a,q+1,r);
		}
	}

	private static int partition(double[] a, int p, int r) {

		double x = a[p];
		int i = p-1 ;
		int j = r+1 ;

		while (true) {
			i++;
			while ( i< r && a[i] < x)
				i++;
			j--;
			while (j>p && a[j] > x)
				j--;

			if (i < j)
				swap(a, i, j);
			else
				return j;
		}
	}

	private static void swap(double[] a, int i, int j) {
		// TODO Auto-generated method stub
		double temp;
		temp = a[i];
		a[i] = a[j];
		a[j] = temp;

	}
	
	//=============================================================

		public static void quickSort(int[] a, int p, int r)
		{
			if(p<r)
			{
				int q=partition(a,p,r);
				quickSort(a,p,q);
				quickSort(a,q+1,r);
			}
		}

		private static int partition(int[] a, int p, int r) {

			double x = a[p];
			int i = p-1 ;
			int j = r+1 ;

			while (true) {
				i++;
				while ( i< r && a[i] < x)
					i++;
				j--;
				while (j>p && a[j] > x)
					j--;

				if (i < j)
					swap(a, i, j);
				else
					return j;
			}
		}

		private static void swap(int[] a, int i, int j) {
			// TODO Auto-generated method stub
			int temp;
			temp = a[i];
			a[i] = a[j];
			a[j] = temp;

		}

}
