/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parseOsm;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.*;

/**
 *
 * @author claudia.plant
 */
public class JSONParser {

    String dir;

    public JSONParser(String dir) {
        this.dir = dir;
        JSONObject data = readFile();
        System.out.println("m");
    }

    private String readFile(String pathname) throws IOException {

    File file = new File(pathname);
    StringBuilder fileContents = new StringBuilder((int)file.length());
    Scanner scanner = new Scanner(file);
    String lineSeparator = System.getProperty("line.separator");

    try {
        while(scanner.hasNextLine()) {        
            fileContents.append(scanner.nextLine() + lineSeparator);
        }
        return fileContents.toString();
    } finally {
        scanner.close();
    }
}
    private JSONObject readFile() {
        // read from the URL
        String data = "";
        try {
            data = readFile(dir);
        } catch (IOException ex) {
            Logger.getLogger(JSONParser.class.getName()).log(Level.SEVERE, null, ex);
        }

        // build a JSON object
        JSONObject obj = new JSONObject(data);
        return obj;
    }
}
