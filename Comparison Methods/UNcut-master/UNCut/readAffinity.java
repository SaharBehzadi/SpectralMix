package graphClustering;

import java.io.*;
import java.util.*;
public class readAffinity {
	ArrayList<double[]> Affinity;
	int numInstance;
	public readAffinity(String filename, int numInstance) throws Exception{
		ArrayList<double[]> Affinity = new ArrayList<double[]>();
		BufferedReader br = new BufferedReader(new FileReader(
				new File(filename)));
		String str[];
		try {
			String line = br.readLine();

			while (line != null) {
				
				str=line.split("\\s+");
				int length = str.length;
				double[] tmpdata = new double[length];
				for (int j=0;j<length;j++){
					tmpdata[j]=Double.parseDouble(str[j]);
				}
				Affinity.add(tmpdata);
				line = br.readLine();
			}
		} finally {
			br.close();
		}
		this.Affinity = Affinity;
		this.numInstance = numInstance;
	}
	
	public ArrayList<double[]> getDataList() {
		return Affinity;
	}
	
	public int[][] getData() {
		int len = Affinity.size();
		int[][] data=new int[len][2];
		int p,q;
		for(int i=0; i<len; i++){
			p=(int)Affinity.get(i)[0];
			q=(int)Affinity.get(i)[1];
			data[i][0]=q;
			data[i][1]=p;
		}

		return data;
	}

}
