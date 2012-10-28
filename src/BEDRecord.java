/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 18.10.12
 * Time: 14:33
 * To change this template use File | Settings | File Templates.
 */

/**
 *
 */
public class BEDRecord {
    private final String contigName;
    private final int stopIndex;
    private final int startIndex;

    /**
     * Immutable record type. Represents a line in BED file.
     * @param contigName Contig name
     * @param startIndex Start index
     * @param stopIndex  Stop index
     */
    public BEDRecord(String contigName, int startIndex, int stopIndex) {
        this.contigName = contigName;
        this.startIndex = startIndex;
        this.stopIndex = stopIndex;
    }

    public int getStartIndex() {
        return startIndex;
    }

    public int getStopIndex() {
        return stopIndex;
    }

    public String getContigName() {
        return contigName;
    }
}
