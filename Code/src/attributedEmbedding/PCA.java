/*
 * PCA.java
 *
 * Created on 25 October 2005, 16:45
 */

package attributedEmbedding;
import Jama.*;


/**
 *
 * @author Administrator
 */
public class PCA {
    DataObject[] data;
    int dim;
    int numObj;
    Matrix reducedVectors;
    
    /** Creates a new instance of PCA */
    public PCA(DataObject[] data) {
        this.data = data;
        dim = data[0].coord.length;
        numObj = data.length;
    }
    
    public DataObject[] pca(int targetDim){
        double[] mean = mean(data);
        double[][] covar = covariance(data, mean);
        EigenvalueDecomposition e = new EigenvalueDecomposition(new Matrix(covar));
        Matrix values = e.getD();
       //values.print(5, 3);
        double[][] rv = new double[targetDim][targetDim];
        for(int i = 0; i < targetDim; i++)
            for(int j = 0; j < targetDim; j++)
                if(i == j)
                    rv[i][j] = values.get(dim-i-1, dim-i-1);
                else
                    rv[i][j] = 0.0;
        Matrix vectors = e.getV();
        //vectors.print(5, 3);
        double[][] rvectors = new double[dim][targetDim];
        for(int i = 0; i < dim; i++)
            for(int j = 0; j < targetDim; j++)
                rvectors[i][j] = vectors.get(i, dim-j-1);
        Matrix reducedVectors = new Matrix(rvectors);
        //reducedVectors.print(5, 3);
        Matrix reducedValues = new Matrix(rv);
        DataObject[] reduced = new DataObject[numObj];
        for(int i = 0; i < numObj; i++){
            double[] coord = new double[targetDim];
            for(int j = 0; j < targetDim; j++){
                for(int l = 0; l < dim; l++)
                    coord[j] += data[i].coord[l] * reducedVectors.get(l, j);
            }
            DataObject red = new DataObject(coord, i);
            red.clusterID = data[i].clusterID;
            red.classID = data[i].classID;
            reduced[i] = red;
        }
        return reduced;
    }
    
    private double[][] covariance(DataObject[] data, double[] mean){
        double[][] cov = new double[dim][dim];
        for(int l = 0; l < data.length; l++)
            for(int i = 0; i < dim; i++)
                for(int j = 0; j < dim; j++)
                    cov[i][j]+= (data[l].coord[j]-mean[j])*(data[l].coord[i]-mean[i]);
        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
                cov[i][j]/= data.length;
        return cov;
    }
    
    public DataObject[] embedPoints(DataObject[] data, int targetDim){
        DataObject[] reduced = new DataObject[data.length];
        for(int i = 0; i < data.length; i++){
            double[] coord = new double[targetDim];
            for(int j = 0; j < targetDim; j++){
                for(int l = 0; l < dim; l++)
                    coord[j] += data[i].coord[l] * reducedVectors.get(l, j);
            }
            DataObject red = new DataObject(coord, i);
            red.clusterID = data[i].clusterID;
            red.classID = data[i].classID;
            reduced[i] = red;
        }
        return reduced;
    }
    
    private double[] mean(DataObject[] members){
        double[] mean = new double[dim];
        for(int i = 0; i < members.length; i++)
            for(int j = 0; j < dim; j++)
                mean[j] += members[i].coord[j]/members.length;
        return mean;
    }
    
}
