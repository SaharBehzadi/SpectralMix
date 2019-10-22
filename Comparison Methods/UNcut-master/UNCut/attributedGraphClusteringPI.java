package graphClustering;

import java.io.*;
import java.util.*;
import smile.math.Math;

public class attributedGraphClusteringPI {

	public static void main(String[] args) throws Exception{
		String dir = "D:/Wei/graph clustering/graph data/genes/poweriter/";
		int numInstance = 2900;
		int numCluster = 15;
//		String dir = args[0];
//		int numInstance = Integer.parseInt(args[1]);

		String fileAffinity=dir+"edgeWeight.txt";
		readAffinity rdf = new readAffinity(fileAffinity,numInstance);
		
//		int[][] Affinity = new int[numInstance][numInstance];
//		Affinity=rdf.getData();
		ArrayList<double[]> edgeWeightList;
		edgeWeightList=rdf.getDataList();
		int num=edgeWeightList.size();
		int[][] edgeWeight = new int[num][3];
		int i,j;
		for (i=0;i<num;i++){
			edgeWeight[i][0]=(int)edgeWeightList.get(i)[0];
			edgeWeight[i][1]=(int)edgeWeightList.get(i)[1];
			edgeWeight[i][2]=(int)edgeWeightList.get(i)[2];
		}
		edgeWeightList=null;
		String fileData=dir+"Data.txt";
		double[][] data;
		readData rdD=new readData(fileData);
		data=rdD.getData();
		
//		String fileLabel=dir+"label.txt";
//		int[] trueLabel;
//		readLabel rdL=new readLabel(fileLabel);
//		trueLabel=rdL.getData();

		long start = System.currentTimeMillis();
		double omega=0.5;

		int[] cluster;
		double phi_star;
		double minobj;

		System.out.println("clustering...");
		
		int[] actuallabel = new int[numInstance];
		
		
		PowerIteration powerIter = new PowerIteration(edgeWeight,data,omega,numCluster,numInstance);
		actuallabel=powerIter.getClusterLabel();

		for (i=0;i<numInstance;i++)
			actuallabel[i]+=1;
		
		String clusterlabel= dir+"clusterlabel.txt";
		PrintWriter clout = new PrintWriter(new FileWriter(clusterlabel));
		for (i = 0; i < numInstance; i++) {
			clout.println(actuallabel[i]);
		}
		clout.close();
		
		long end = System.currentTimeMillis();
		long difference = end - start;
//		String filename=dir+ String.valueOf(omega);
//		File f = new File(filename);
//        f.mkdir();
//		String clusterlabel= dir+String.valueOf(omega) + "/clusterlabel.txt";
//		PrintWriter clout = new PrintWriter(new FileWriter(clusterlabel));
//		for (i = 0; i < numInstance; i++) {
//			clout.println(actuallabel[i]);
//		}
//		clout.close();
		
		System.out.println("eclipsed time: "+ difference/1000.0 +'s');

		String results= dir + "/results.txt";
		PrintWriter rout = new PrintWriter(new FileWriter(results));
		rout.print("eclipsed time: ");
		rout.print(difference/1000.0);
		rout.println("sec");
		rout.close();

		System.out.println("Done!");

	}

}

