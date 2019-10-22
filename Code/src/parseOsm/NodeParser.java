/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parseOsm;

//import br.zuq.osm.parser.model.OSMNode;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

/**
 *
 * @author zuq
 */
public class NodeParser {

    public static boolean isNode(Node node) {
        return (node.getNodeName().equals("node"));
    }

    public static OSMNode parseNode(Node node) {
        NamedNodeMap atts = node.getAttributes();

        String id = atts.getNamedItem("id").getNodeValue();

        OSMNode osmNode = new OSMNode(id,
                atts.getNamedItem("visible").getNodeValue(),
                atts.getNamedItem("timestamp").getNodeValue(),
                atts.getNamedItem("version").getNodeValue(),
                atts.getNamedItem("changeset").getNodeValue(),
                atts.getNamedItem("user").getNodeValue(),
                atts.getNamedItem("uid").getNodeValue(),
                atts.getNamedItem("lat").getNodeValue(),
                atts.getNamedItem("lon").getNodeValue(),
                OSMParser.parseTags(node.getChildNodes()));
       
        return osmNode;
    }
}
