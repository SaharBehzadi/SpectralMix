package parseOsm;

/**
 *
 * @author Willy Tiengo
 */
public class AbstractNode {

    public String id;
    public String visible;
    public String timestamp;
    public String version;
    public String changeset;
    public String user;
    public String uid;

    public AbstractNode(String id, String visible, String timestamp,
            String version, String changeset, String user, String uid) {

        this.id = id;
        this.visible = visible;
        this.timestamp = timestamp;
        this.version = version;
        this.changeset = changeset;
        this.user = user;
        this.uid = uid;

    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }

        final OSMNode other = (OSMNode) obj;

        return (id != null ? id.equals(other.id) : false);

    }

    @Override
    public int hashCode() {
        return (id != null ? id.hashCode() : super.hashCode());
    }
}

