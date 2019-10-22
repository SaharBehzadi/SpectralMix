package graphClustering;

import smile.clustering.KMeans;

import java.io.*;
import java.util.*;

public class PowerIteration {

	int[] clusterLabel;


	public PowerIteration(int[][] A, double[][] data, double omega, int numCluster, int n) throws Exception{

		int dim=data[0].length;
		int i,j;
		int candik=10*numCluster;

		double[] D = new double[n];
		int[] p = new int[n];
		int[] index = new int[n+1];
		index[0]=0;
		index[n]=A.length;

		for (i = 1; i < n; i++) {
			for (j = 0; j < A.length-1; j++) {
				if ((A[j][0]==i-1)&&(A[j+1][0]==i)){
					index[i]=j+1;
					break;
				}
			}
		}

		for (i=0; i<n; i++)
			D[i]=index[i+1]-index[i];


		double[][] W = new double[A.length][3];
		for (i = 0; i < A.length; i++) {
			for (j = 0; j < 3; j++) {
				W[i][j] = A[i][j];
			}
		}

		for (i=0; i<n;i++){
			for (j=index[i]; j<index[i+1];j++)
				W[j][2]=1.0*A[j][2]/D[i];
		}

		double[][] L = new double[A.length][3];
		for (i=0; i<A.length; i++){
			L[i][0]=W[i][0];
			L[i][1]=W[i][1];
			L[i][2]=-W[i][2];
		}

		double[][] Y=new double[candik][n];
		double[] lamda = new double[candik];
		double[] u = new double[n];
		double tmp=0;

		for (i=0; i<candik; i++){

			PI powerIter = new PI(W,n,index);
			Y[i]=powerIter.getPseudoEig();
			for (int i1=0; i1<n; i1++){
				tmp=0;
				tmp+=Y[i][i1];
				for (j=index[i1]; j<index[i1+1];j++){
					int col=(int)L[j][1];
					tmp += L[j][2]*Y[i][col];
				}
				u[i1]=tmp;
			}
			tmp=0;
			for (int i1=0; i1<n;i1++)
				tmp += u[i1]*Y[i][i1];

			lamda[i]=tmp;


		}

		long start = System.currentTimeMillis();
		double[][] eigenvector = new double[n][2];
		int[] label = new int[n];
		int cnt1=0;
		int cnt2=0;
		ArrayList<Integer> index1 = new ArrayList<Integer>();
		ArrayList<Integer> index2 = new ArrayList<Integer>();
		double[] unicomp = new double[candik];
		double[][] subdata;
		for (i=0;i<candik;i++)
			unicomp[i]=1000;

		for (i=0;i<candik;i++){
			for (j=0;j<n;j++){
				eigenvector[j][0]=Y[i][j];
				eigenvector[j][1]=0;
			}

			KMeans kmeans = new KMeans(eigenvector, 2);
			label=kmeans.getClusterLabel();
			cnt1=0;
			cnt2=0;
			for (j=0;j<n;j++){
				if (label[j]==0){
					cnt1++;
					index1.add(j);
				}
				else{
					cnt2++;
					index2.add(j);
				}

			}
			if (cnt1<cnt2){
				if (cnt1==0)
					unicomp[i]=1000;
				else{
					subdata = new double[cnt1][dim];
					for (int k=0;k<cnt1;k++){
						for (j=0;j<dim;j++)
							subdata[k][j]=data[index1.get(k)][j];
					}
					unimodCompact compactness = new unimodCompact(subdata);
					unicomp[i]=compactness.getUniComp();
				}
			}
			if (cnt1>cnt2){
				if (cnt2==0)
					unicomp[i]=1000;
				else{
					subdata = new double[cnt2][dim];
					for (int k=0;k<cnt2;k++){
						for (j=0;j<dim;j++)
							subdata[k][j]=data[index2.get(k)][j];
					}
					unimodCompact compactness = new unimodCompact(subdata);
					unicomp[i]=compactness.getUniComp();
				}
			}
		}
		long end = System.currentTimeMillis();
		long difference = end - start;
		System.out.println("eclipsed time for the dip test: "+ difference/1000.0 +'s');

		double[][] phi = new double[2][candik];
		for (i=0;i<candik;i++){
			phi[0][i]=(1-omega)*lamda[i]+omega*unicomp[i];
			phi[1][i]=i;
		}




		myQuickSort.quickSort(phi, 0, candik-1);

		//		for (int k1=0;k1<numEigenVec;k1++){
		//			System.out.print(lamda[k1][k1]);
		//			System.out.print(' ');
		//		}
		//		System.out.println();
		//		for (int k1=0;k1<numEigenVec;k1++){
		//			System.out.print(unicomp[k1]);
		//			System.out.print(' ');
		//		}
		//		System.out.println();


		double[][] resultEigenvector = new double[n][numCluster];

		for (i=0;i<numCluster;i++){
			for (j=0;j<n;j++){
				resultEigenvector[j][i]=Y[(int)phi[1][i]][j];
			}
		}

		int[] clusterLabel = new int[n]; 
		KMeans kmeans = new KMeans(resultEigenvector, numCluster);
		clusterLabel=kmeans.getClusterLabel();
		this.clusterLabel=clusterLabel;
	}

	public int[] getClusterLabel() {
		return clusterLabel;
	}

}
