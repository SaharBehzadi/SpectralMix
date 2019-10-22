package graphClustering;

import java.util.*;

public class unimodCompact{
	double uniComp;
	int numMulticluster;

	public unimodCompact(double[][] data) throws Exception{
		int cnt=0,i,j,sum=0;
		double significance=0.05;
		double uniComp;
		int d=data[0].length,n=data.length;

		int numMulticluster=0;
		
		double[] tmpdata=new double[n];
		ArrayList<Double> dip=new ArrayList<Double>();

		for (i=0;i<d;i++){
			for (j=0;j<n;j++)
				tmpdata[j]=data[j][i];
			myQuickSort.quickSort(tmpdata,0, n-1);
			
			hartigansDipSignifTest Diptest = new hartigansDipSignifTest(tmpdata);

			if (Diptest.getpValue()>=significance){
				dip.add(Diptest.getDip());
			}
			else
				cnt++;
		}

		if (cnt==d){
			uniComp=2*Math.log(d)/Math.log(2);
//			uniComp=2*d;
			numMulticluster=1;
		}
		else{
			Iterator<Double> iter = dip.iterator();
			while(iter.hasNext()){
				sum+=iter.next();
			}
			uniComp=Math.log(1.0*d/dip.size())/Math.log(2)+1.0*sum/dip.size();
//			uniComp=1.0*d/dip.size()+1.0*sum/dip.size();
		}
		
		this.uniComp=uniComp;
		this.numMulticluster=numMulticluster;
	}
	
	public double getUniComp() {
		return uniComp;
	}
	public double getMultimodal() {
		return numMulticluster;
	}

}
