package nature;
import java.io.File;
import java.io.FileWriter;
import weka.core.Instance;

import weka.core.Instances;
import weka.core.Utils;

public class ArffFileWriter {

	public ArffFileWriter() {
		
	}

	public boolean saveFile(String path, Instances instances) {
		
		System.out.println("saving dataset: " + path);
		File file = new File(path);
		
		try {
			FileWriter fileWriter = new FileWriter(file);
			fileWriter.write("@relation test\n\n");
			
			for (int i = 0; i < instances.numAttributes(); i++) {
				fileWriter.write(instances.attribute(i).toString() + "\n");
			}
			fileWriter.write("\n");
			fileWriter.write("@data\n");
                        
			for (int i = 0; i < instances.numInstances(); i++) {
				fileWriter.write(instances.instance(i).toString() + "\n");
                               
			}
			fileWriter.close();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return true;
	}
        
        //avoid the usage of weka.core.instance.toString() since this causes problems in some cases
        public boolean saveFileWithoutInstanceToString(String path, Instances instances){
            System.out.println("saving dataset: " + path);
		File file = new File(path);
		
		try {
			FileWriter fileWriter = new FileWriter(file);
			fileWriter.write("@relation test\n\n");
			
			for (int i = 0; i < instances.numAttributes(); i++) {
				fileWriter.write(instances.attribute(i).toString() + "\n");
			}
			fileWriter.write("\n");
			fileWriter.write("@data\n");
                     
			for (int i = 0; i < instances.numInstances(); i++) {
                            StringBuffer text = new StringBuffer();
                            for(int j = 0; j < instances.numAttributes(); j++){
                                   if (j > 0) text.append(",");
                                 text.append(Utils.doubleToString(instances.instance(i).value(j),6));
                            }
				fileWriter.write(text + "\n");
                               
			}
			fileWriter.close();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return true;
            
        }

}

