/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parseOsm;

//import br.zuq.osm.parser.model.Member;
//import br.zuq.osm.parser.model.Relation;
import java.util.ArrayList;
import java.util.List;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
 *
 * @author zuq
 */
public class RelationParser {

    public static boolean isRelation(Node node) {
        return node.getNodeName().equals("relation");
    }

    public static Relation parseRelation(Node node) {
        NamedNodeMap atts = node.getAttributes();

        String id = atts.getNamedItem("id").getNodeValue();

        return new Relation(id,
                atts.getNamedItem("visible").getNodeValue(),
                atts.getNamedItem("timestamp").getNodeValue(),
                atts.getNamedItem("version").getNodeValue(),
                atts.getNamedItem("changeset").getNodeValue(),
                atts.getNamedItem("user").getNodeValue(),
                atts.getNamedItem("uid").getNodeValue(),
                getMembers(node.getChildNodes()),
                OSMParser.parseTags(node.getChildNodes()));
    }

    private static List<Member> getMembers(NodeList children) {
        List<Member> result;
        Node node;
        NamedNodeMap map;

        result = new ArrayList<Member>();

        for (int i = 0; i < children.getLength(); i++) {
            node = children.item(i);
            map = node.getAttributes();

            if (node.getNodeName().equals("member")) {
                result.add(new Member(
                        map.getNamedItem("type").getNodeValue(),
                        map.getNamedItem("ref").getNodeValue(),
                        map.getNamedItem("role").getNodeValue()));
            }
        }

        return result;
    }
}

