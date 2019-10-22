package graphClustering;

import java.util.ArrayList;

public class realdataQuality {
	public static void main (String[] args) throws Exception{
		String dir = "C:/Users/ye/Desktop/graph clustering/graph data/dfb/";
		int numCluster=6;
		
		String fileclusterlabel=dir+"clusterlabel.txt";
		readLabel rd = new readLabel(fileclusterlabel);
		int[] uncutlabel;
		uncutlabel=rd.getData();
		
		fileclusterlabel=dir+"sscglabel.txt";
		rd = new readLabel(fileclusterlabel);
		int[] sscglabel;
		sscglabel=rd.getData();
		
		fileclusterlabel=dir+"SAcluster/asignment.txt";
		rd = new readLabel(fileclusterlabel);
		int[] salabel;
		salabel=rd.getData();
		
		int numInstance=uncutlabel.length;
		String fileAffinity=dir+"edgeWeight.txt";
		readAffinity rdf = new readAffinity(fileAffinity,numInstance);
		int[][] Affinity = new int[numInstance][numInstance];
		Affinity=rdf.getData();
		String fileData=dir+"Data.txt";
		double[][] data;
		readData rdD=new readData(fileData);
		data=rdD.getData();
		
		double[] uncutResult = new double[2];
		double[] sscgResult = new double[2];
		double[] saResult = new double[2];
		uncutResult=uncutComputation(data,Affinity,uncutlabel,numCluster);
		sscgResult=uncutComputation(data,Affinity,sscglabel,numCluster);
		saResult=uncutComputation(data,Affinity,salabel,numCluster);
		
		System.out.println("uncut:");
		System.out.print(uncutResult[0]);
		System.out.print('\t');
		System.out.println(uncutResult[1]);
		
		System.out.println("sscg:");
		System.out.print(sscgResult[0]);
		System.out.print('\t');
		System.out.println(sscgResult[1]);
		
		System.out.println("sacluster:");
		System.out.print(saResult[0]);
		System.out.print('\t');
		System.out.println(saResult[1]);
	}
	
	public static double[] uncutComputation(double[][] data, int[][] Affinity, int[] label, int numCluster) throws Exception{
		int dim=data[0].length;
		int n=data.length;
		ArrayList<Integer> index = new ArrayList<Integer>();
		double unicomp =0;
		int numMultimodalCluster=0;
		int[] one = new int[n];
		for (int i=0;i<n;i++)
			one[i]=1;
		
		double ncut=0;
		
		for (int k=1;k<=numCluster;k++){
			//compute unimodal compatness
			int cnt=0;
			for (int i=0;i<label.length;i++){
				if (label[i]==k){
					cnt++;
					index.add(i);
				}
			}
			int[] B = new int[n];
			double[][] subdata = new double[cnt][dim];
			
			for (int i=0;i<cnt;i++){
				int c=index.get(i);
				B[c]=1;
				for (int j=0;j<dim;j++)
					subdata[i][j]=data[c][j];
			}
			
			unimodCompact compactness = new unimodCompact(subdata);
			unicomp += compactness.getUniComp();
			if (compactness.getMultimodal()==1)
				numMultimodalCluster++;
			
			//compute normalized cut
			int[] sub = new int[n];
			for (int i=0;i<n;i++)
				sub[i]=one[i]-B[i];
			
			double[] multi = new double[n];
			double tmp=0;
			for (int i=0;i<n;i++){
				tmp=0;
				for (int j=0;j<n;j++){
					tmp += B[j]*Affinity[j][i];
				}
				multi[i]=tmp;
			}
			double numerator=0, denominator=0;
			for (int i=0;i<n;i++){
				numerator += multi[i]*sub[i];
				denominator += multi[i]*one[i];
			}
			
			ncut += numerator/denominator;		
		}
//		System.out.print("ncut: ");
//		System.out.println(ncut);
//		System.out.print("unicomp: ");
//		System.out.println(unicomp);
		
		double[] result = new double[2];
		result[0]=ncut+unicomp;
		result[1]=numMultimodalCluster;
		
		return result;
		
	}

}
