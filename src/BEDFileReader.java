import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 18.10.12
 * Time: 14:27
 * To change this template use File | Settings | File Templates.
 */
public class BEDFileReader {
    BEDFileReader(File in)
    {
        Logger logger = Logger.getLogger(this.getClass());
        try{
            this.input = new BufferedReader(new FileReader(in));
            this.buffer = this.input.readLine();
        } catch (Exception e){
            logger.error("Error opening BED file: " + e.toString());
        }
    }

    public BEDRecord getBEDRecord()
    {
        Logger logger = Logger.getLogger(this.getClass());
        try {
            String[] data = buffer.split("\t");
            buffer = input.readLine();
            return new BEDRecord(data[0], Integer.parseInt(data[1]), Integer.parseInt(data[2]));
        } catch (IOException e) {
            logger.error("Error reading BED file: " + e.toString());
            return new BEDRecord("null", 0, 0);
        }
    }
    public boolean hasNext()
    {
        return (null != buffer);
    }

    public void close(){
        Logger logger = Logger.getLogger(this.getClass());
        try {
            input.close();
        } catch (IOException e) {
            logger.error("Error closing BED file: " + e.toString());
        }
    }

    private BufferedReader input = null;
    private String buffer = null;
}
