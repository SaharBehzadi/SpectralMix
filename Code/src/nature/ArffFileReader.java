/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nature;
import java.io.File;
import java.io.FileReader;

import weka.core.Instances;

/**
 *
 * @author claudia.plant
 */
public class ArffFileReader {
  

	public Instances readFile(String filePath) {
		
		//System.out.println("loading dataset...");
		File file = new File(filePath);
		
		try {
			FileReader fileReader = new FileReader(file);
			Instances instances = new Instances(fileReader);
			instances.setClassIndex(instances.numAttributes() - 1);
            //System.out.println(instances.numInstances() + " instances with "
            //		+ instances.numAttributes() + " attributes loaded.");
            fileReader.close();
			return instances;
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return null;
	}

}


    

