/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;

import com.jmatio.io.MatFileReader;
import com.jmatio.io.MatFileWriter;
import com.jmatio.types.MLArray;
import com.jmatio.types.MLDouble;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.InstanceComparator;
import weka.core.Instances;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.NumericToNominal;
import weka.filters.unsupervised.attribute.Standardize;

/**
 *
 * @author claudia
 */
public class DataTransformer {

    Instances data;

    public DataTransformer(Instances data) {
        this.data = data;
    }

    public DataTransformer() {
    }

    public Instances zScale(Instances dd) {
        try {
            Standardize filter = new Standardize();
            filter.setIgnoreClass(true);
            filter.setInputFormat(dd);
            dd = Filter.useFilter(dd, filter);
        } catch (Exception ex) {
            Logger.getLogger(DataTransformer.class.getName()).log(Level.SEVERE, null, ex);
        }
        return dd;
    }

    public Instances merge(Instances a, Instances b) {
        Instances m = new Instances(a);
        for (int i = 0; i < b.numInstances(); i++) {
            m.add(b.instance(i));
        }
        return m;
    }

    public double[][] extract(Instances d) {
        double[][] dd = new double[d.numInstances()][d.numAttributes()];
        for (int i = 0; i < d.numInstances(); i++) {
            for (int j = 0; j < d.numAttributes(); j++) {
                dd[i][j] = d.instance(i).value(j);
            }
        }
        return dd;
    }

    public DataObject[] createD(double[][] coord) {
        int dim = coord[0].length;
        int numObj = coord.length;
        DataObject[] res = new DataObject[numObj];
        for (int i = 0; i < res.length; i++) {
            res[i] = new DataObject(coord[i], i);
        }
        return res;
    }

    public DataObject[] createD1d(double[] coord) {
        int numObj = coord.length;
        DataObject[] res = new DataObject[numObj];
        for (int i = 0; i < res.length; i++) {
            double[] c = new double[1];
            c[0] = coord[i];
            res[i] = new DataObject(c, i);
        }
        return res;
    }

    public Instances create(double[] coord) {
        int dim = 1;
        int numObj = coord.length;
        FastVector fvWekaAttributes = new FastVector(dim);

        fvWekaAttributes.addElement(new Attribute("bla"));

        Instances t = new Instances("test", fvWekaAttributes, numObj);
        for (int i = 0; i < numObj; i++) {
            Instance inst = new Instance(dim);
            for (int j = 0; j < dim; j++) {
                inst.setValue(j, coord[i]);
            }
            t.add(inst);
        }
        return t;
    }

    public Instances removeMissing(Instances inst) {
        Instances dataWMissing = new Instances(inst);
        dataWMissing.delete();
        for (int i = 0; i < inst.numInstances(); i++) {
            boolean missing = false;
            for (int j = 0; j < inst.numAttributes(); j++) {
                if (inst.instance(i).isMissing(j)) {
                    missing = true;
                }

            }
            if (!missing) {
                dataWMissing.add(inst.instance(i));
            }
        }
        return dataWMissing;
    }

    public void writeArff(double[][] coord, int[] clusterId, int[] classId, String path) {
        try {
            int dim = coord[0].length;
            int numObj = coord.length;
            FastVector fvWekaAttributes = new FastVector(dim + 2);
            for (int i = 0; i < dim; i++) {
                fvWekaAttributes.addElement(new Attribute("bla" + Integer.toString(i)));
            }
            fvWekaAttributes.addElement(new Attribute("clusterId"));
            fvWekaAttributes.addElement(new Attribute("classId"));
            Instances t = new Instances("test", fvWekaAttributes, numObj);
            for (int i = 0; i < numObj; i++) {
                Instance inst = new Instance(dim);
                for (int j = 0; j < dim; j++) {
                    inst.setValue(j, coord[i][j]);
                }
                inst.setValue(dim, clusterId[i]);
                inst.setValue(dim + 1, classId[i]);
                t.add(inst);
            }
            NumericToNominal n = new NumericToNominal();
            int[] atts = new int[2];
            atts[0] = dim;
            atts[1] = dim + 1;
            n.setAttributeIndicesArray(atts);
            n.setInputFormat(t); // inform filter about dataset **AFTER** setting options
            t = Filter.useFilter(t, n);
            System.out.println("saving dataset: " + path);
            File file = new File(path);
            try {
                FileWriter fileWriter = new FileWriter(file);
                fileWriter.write("@relation test\n\n");
                for (int i = 0; i < t.numAttributes(); i++) {
                    fileWriter.write(t.attribute(i).toString() + "\n");
                }
                fileWriter.write("\n");
                fileWriter.write("@data\n");
                for (int i = 0; i < t.numInstances(); i++) {
                    fileWriter.write(t.instance(i).toString() + "\n");
                }
                fileWriter.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } catch (Exception ex) {
            Logger.getLogger(DataTransformer.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public Instances create(double[][] coord) {
        int dim = coord[0].length;
        int numObj = coord.length;
        FastVector fvWekaAttributes = new FastVector(dim);
        for (int i = 0; i < dim; i++) {
            fvWekaAttributes.addElement(new Attribute("bla" + Integer.toString(i)));
        }
        Instances t = new Instances("test", fvWekaAttributes, numObj);
        for (int i = 0; i < numObj; i++) {
            Instance inst = new Instance(dim);
            for (int j = 0; j < dim; j++) {
                inst.setValue(j, coord[i][j]);
            }
            t.add(inst);
        }
        return t;
    }

    public double[][] readMatlabMatrix(String dir, String variableName) {
        double[][] dataM = new double[1][1];

        MatFileReader mfr = null;
        try {
            mfr = new MatFileReader(dir);
        } catch (IOException e) {
        }

        if (mfr != null) {
            dataM = ((MLDouble) mfr.getMLArray(variableName)).getArray();
        }
        return dataM;
    }

    public void saveAsMatlab(double[][] d, String vName, String fName) {
        MLDouble v = new MLDouble(vName, d);
        ArrayList ll = new ArrayList();
        ll.add(v);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        System.out.println("Matlab matrix saved.");
    }

    public void saveAsMatlab2(double[][] d, double[][] e, String vName_d, String vName_e, String fName) {
        MLDouble v = new MLDouble(vName_d, d);
        MLDouble v1 = new MLDouble(vName_e, e);
        ArrayList ll = new ArrayList();
        ll.add(v);
        ll.add(v1);
        MatFileWriter mw = new MatFileWriter();
        try {
            mw.write(fName, ll);
        } catch (IOException ex) {
        }
        System.out.println("Matlab matrix saved.");
    }

    //29.12.2011
    public void saveSuperIndicatorMatrix(Instances data) {
        int totalNumValues = 0;
        for (int i = 0; i < data.numAttributes(); i++) {
            if (data.attribute(i).isNominal()) {
                totalNumValues += data.attribute(i).numValues();
            }
        }
        double[][] ind = new double[data.numInstances()][totalNumValues];
        for (int i = 0; i < data.numInstances(); i++) {
            int offSet = 0;
            for (int j = 0; j < data.numAttributes(); j++) {
                ind[i][(int) data.instance(i).value(j) + offSet] = 1.0;
                offSet += data.attribute(j).numValues();
            }
        }
        saveAsMatlab(ind, "superInd", "superInd.mat");
    }

    public Instances getCluster(Instances d, double[] id, int label) {
        Instances res = new Instances(d);
        res.delete();
        for (int i = 0; i < d.numInstances(); i++) {
            if (id[i] == label) {
                res.add(new Instance(d.instance(i)));
            }
        }

        return res;

    }

    public void extractAttInformation(Instances d, int index, int numClusters) {
        for (int i = 0; i < numClusters; i++) {
        }
    }

    public Instances removeDuplicates(Instances d) {
        Instances toSort = new Instances(d);
        Instance[] toSortA = new Instance[toSort.numInstances()];
        for (int i = 0; i < toSort.numInstances(); i++) {
            toSortA[i] = toSort.instance(i);
        }
        InstanceComparator ic = new InstanceComparator();
        ic.setIncludeClass(false);
        Arrays.sort(toSortA, ic);
        Instances res = new Instances(d);
        res.delete();
        for (int i = 0; i < toSortA.length; i++) {
            res.add(toSortA[i]);
        }

        boolean finished = false;
        int self = 0;
        int compare = 1;
        while (!finished) {
            if (ic.compare(res.instance(self), res.instance(compare)) == 0) {
                res.delete(compare);
                self = 0;
                compare = 1;

            } else {
                self++;
                compare++;
            }
            if (self == res.numInstances() - 1) {
                finished = true;
            }

        }
        return res;

    }

    public Instances[] separateClusters(Instances d, double[] id, int num) {
        Instances[] result = new Instances[num];
        for (int i = 0; i < result.length; i++) {
            result[i] = new Instances(d);
            result[i].delete();
        }
        for (int i = 0; i < d.numInstances(); i++) {
            result[(int) id[i] - 1].add(new Instance(d.instance(i)));
        }
        return result;
    }
}
