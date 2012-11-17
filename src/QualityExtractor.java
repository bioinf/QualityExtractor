import com.sun.deploy.util.ArrayUtil;
import net.sf.samtools.*;
import org.apache.log4j.Logger;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.*;

public class QualityExtractor {

    /**
     *
     * @param args arg[0] is BAM file, args[1] is BED file
     */
    public static void main(String[] args) {

        Logger logger = Logger.getLogger(QualityExtractor.class);

        if(args.length != 2) {
            logger.error("Not enough arguments to launch. Quit.");
            return;
        }

        logger.info("Start processing...");
        createBamIndex(new File(args[0]), new File(args[0] + ".bai"));
        calcMeanQuality(extractReadsFromBAM(new File(args[0]), new File(args[1])));
        logger.info("Processing finished.");
    }

    /**
     * Create a BAM index file. If output file with this name already exists, no index file will be created
     * @param inputSamOrBamFile   Input BAM file name
     * @param outputBamIndexFile    Output BAM index file name
     */
    public static void createBamIndex(final File inputSamOrBamFile, final File outputBamIndexFile)
    {
        Logger logger = Logger.getLogger(QualityExtractor.class);
        try{
            logger.info("BAM indexing started");
            if(outputBamIndexFile.exists())
            {
                logger.info("BAM index file " + outputBamIndexFile.getName() + " already exists. New index will not be created.");
                return;
            }
            final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
            BAMIndexer indexer = new BAMIndexer(outputBamIndexFile, inputSam.getFileHeader());
            inputSam.enableFileSource(true);
            int totalRecords = 0;

            for(SAMRecord rec : inputSam)
            {
                if (++totalRecords % 10000 == 0)
                    logger.info(totalRecords + " reads processed ...");

                indexer.processAlignment(rec);
            }
            indexer.finish();
            inputSam.close();

            logger.info("BAM indexing finished");
        }catch (Exception e){
            logger.error("Error while data processing: " + e.toString());
        }
    }

    /**
     * Extract all reads, that fit the boundaries given by BED file
     * @param inputSamOrBamFile
     * @param inputBedFile
     * @return Map with BEDRecords and associated non-empty list of SAMRecords. If BEDRecord has no associated
     * SAMRecords, it will not be in the result map
     */
    public static  Map<BEDRecord, List<SAMRecord>> extractReadsFromBAM(final File inputSamOrBamFile, final File inputBedFile) {
        Logger logger = Logger.getLogger(QualityExtractor.class);
        Map<BEDRecord, List<SAMRecord>> resDict = new HashMap<BEDRecord, java.util.List<SAMRecord>>();
        try{
            final BEDFileReader inputBed = new BEDFileReader(inputBedFile);
            final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile, new File(inputSamOrBamFile.getName() + ".bai"));
            if(!inputSam.hasIndex())
                throw new Exception("Index file is invalid");
            while (inputBed.hasNext())
            {
                BEDRecord record = inputBed.getBEDRecord();
                logger.info("BEDRecord: " + record.getContigName()+ "\t" + record.getStartIndex() + '\t' + record.getStopIndex());
                SAMRecordIterator it = inputSam.query(record.getContigName(), record.getStartIndex(), record.getStopIndex(), false);
                List<SAMRecord> list = new ArrayList<SAMRecord>();
                int totalRecords = 0;
                while (it.hasNext())
                {
                    SAMRecord rec = it.next();
                    list.add(rec);
                    totalRecords++;
                    //logger.info("SAMRecord: " + rec.getReadString());
                }
                it.close();
                logger.info("Total SAMRecords for current BED: " + totalRecords);

                if(list.size() > 0)
                    resDict.put(record, list);
            }
            inputBed.close();
            inputSam.close();

        }catch (Exception e){
            logger.error("Error iterating through reads: " + e.toString());
        }
        return  resDict;
    }

    /**
     * Parse the CIGAR string for the record, and return aligned quality and Read arrays
     * @param record
     * @return
     */
    public static void getAlignedQualities(SAMRecord record, byte[] outputQuality, byte[] outputRead)
    {
        Logger logger = Logger.getLogger(QualityExtractor.class);
        logger.debug("Read Name: " + record.getReadName());
        logger.debug("Read string: " + record.getReadString());
        logger.debug("Quality string: " + record.getBaseQualityString());
        logger.debug("Cigar string: " + record.getCigarString());
        logger.debug("Quality array size: " + record.getReadBases().length + "; It starts from " +
                record.getAlignmentStart() + " and ends at " + record.getAlignmentEnd() + "; Difference of alignment " +
                "positions gives us length " + (record.getAlignmentEnd() - record.getAlignmentStart() + 1));

        ArrayList<Byte> resultRead = new ArrayList<Byte>();
        ArrayList<Byte> resultQuality = new ArrayList<Byte>();
        List<CigarElement> elements = record.getCigar().getCigarElements();

        //convert quality and read arrays to stack for convinience
        ArrayList<Byte> read = new ArrayList<Byte>();
        for(byte b : record.getReadBases())
            read.add(new Byte(b));
        Collections.reverse(read);
        Stack<Byte> reads = new Stack<Byte>();
        reads.addAll(read);

        ArrayList<Byte> qual = new ArrayList<Byte>();
        for(byte b : record.getBaseQualities())
            qual.add(new Byte(b));
        Collections.reverse(qual);
        Stack<Byte> qualities = new Stack<Byte>();
        qualities.addAll(qual);

        //TODO think of how to use AlignmentBlock instead
        for(CigarElement elem : elements) {
            switch (elem.getOperator()) {
                case M:
                    for(int i = 0; i < elem.getLength(); ++i){
                        resultRead.add(reads.pop());
                        resultQuality.add(qualities.pop());
                    }
                    break;
                case D:
                    for(int i = 0; i < elem.getLength(); ++i){
                        resultRead.add(new Byte((byte)'-'));
                        resultQuality.add(new Byte((byte)'*'));
                    }
                    break;
                case S:
                case I:
                    for(int i = 0; i < elem.getLength(); ++i){
                        reads.pop();
                        qualities.pop();
                    }
                    break;
                case N:
                case H:
                case P:
                case EQ:
                case X:
                    logger.error("Unsupported CIGAR operator for read " + record.getReadName());
            }
        }

        for(int i = 0; i < resultRead.size(); ++i){
            outputRead[i] = resultRead.get(i);
            outputQuality[i] = resultQuality.get(i);
        }

        //deal with DEL - replace quality with the mean value from the neighbours
        int i = 0;
        while(i < outputQuality.length)
        {
            if(outputQuality[i] != '*')
            {
                i++;
                continue;
            }
            if(outputQuality.length == i)
                break;
            int left_boundary = i;
            while(i < outputQuality.length && outputQuality[i] == '*')
                i++;
            int right_boundary = i - 1;

            byte mean_quality;
            if(0 == left_boundary && outputQuality.length - 1 == right_boundary)
                mean_quality = 0;
            else if(0 == left_boundary)
                mean_quality = outputQuality[right_boundary + 1];
            else if(outputQuality.length - 1 == right_boundary)
                mean_quality = outputQuality[left_boundary - 1];
            else //normal case - we have neighbours on the left and on the right
                mean_quality = (byte)((outputQuality[left_boundary - 1] + outputQuality[right_boundary + 1]) / 2);

            //fill DEL regions with mean values from the neighbours
            for(int j = 0; j < right_boundary - left_boundary + 1; ++j)
                outputQuality[left_boundary + j] = mean_quality;
        }
        logger.debug("Aligned read and quality strings of length " + outputRead.length + ": " + new String(outputRead).toString());
    }

    /**
     * Calculates mean quality for particular BEDRecord's associated reads in SAMRecord
     * @param dict BEDRecords with associated SAMRecords
     */
    public static void calcMeanQuality(Map<BEDRecord, List<SAMRecord>> dict) {
        Logger logger = Logger.getLogger(QualityExtractor.class);

        try{
            Set s = dict.entrySet();
            Iterator it = s.iterator();

            //Iterate through map records
            while(it.hasNext())
            {
                Map.Entry m =(Map.Entry)it.next();
                BEDRecord bedRecord = (BEDRecord) m.getKey();
                List<SAMRecord> records = (List<SAMRecord>) m.getValue();
                logger.info("Parsing reference sequence: " + bedRecord.getContigName() + " " + bedRecord.getStartIndex() +
                        " " + bedRecord.getStopIndex() + " with " + records.size() + " sequences");

                //three-dimensional array for calculating mean value for all reads int this position
                //0-th dimension - nucleotide (A, C, G, T, DEL). We calculate quality separately for each nucleotide
                //1-st dimension - sum of qualities of all reads for i-th position in the aligned sequence,
                //2-nd dimension - number of reads that were used to calculate the value in the first row
                //So, totalQuality[i][0] / totalQuality[i][1] is the mean value for quality for this position
                int [][][] totalQuality = new int[5][bedRecord.getStopIndex() - bedRecord.getStartIndex() + 1][2];

                logger.info("Start calculating the mean value for quality...");
                for(SAMRecord samRecord : records){
                    int samStartPos = samRecord.getAlignmentStart(), samStopPos = samRecord.getAlignmentEnd(),
                            bedStartPos = bedRecord.getStartIndex(), bedStopPos = bedRecord.getStopIndex();

                    if(0 == samStartPos) {
                        continue;     //unaligned read. Skip
                    }

                    byte[] quality = new byte[samStopPos - samStartPos + 1];
                    byte[] read = new byte[samStopPos - samStartPos + 1];
                    QualityExtractor.getAlignedQualities(samRecord, quality, read);

                    //find indicies to correctly copy overlapped record
                    int i = 0, // start index in SAMRecord
                            j = 0, // start index in BEDRecord
                            copySize = 0; //num of bytes to copy from SAMRecord
                    //SAMRecord can overlap BEDRecord, so we need to be careful
                    if(samStartPos >= bedStartPos && samStopPos <= bedStopPos){
                        //SAMRecord is included in BED boundaries
                        i = 0;
                        j = samStartPos - bedStartPos;
                        copySize = samStopPos - samStartPos + 1;
                    } else if(samStartPos < bedStartPos && samStopPos <= bedStopPos){
                        //SAMRecord overlaps the left boundary
                        i = bedStartPos - samStartPos;
                        j = 0;
                        copySize = samStopPos - bedStartPos + 1;
                    } else if(samStartPos >= bedStartPos && samStopPos > bedStopPos){
                        //SAMRecord overlaps the right boundary
                        i = 0;
                        j = samStartPos - bedStartPos;
                        copySize = bedStopPos - samStartPos + 1;
                    } else if(samStartPos < bedStartPos && samStopPos > bedStopPos){
                        //SAMRecord overlaps both boundaries
                        i = bedStartPos - samStartPos;
                        j = 0;
                        copySize = bedStopPos - bedStartPos + 1;
                    }

                    //copy SAMRecord to BEDRecord boundaries
                    copySize += i; //if i != 0, because we need to iterate copySize times
                    while(i < copySize)
                    {
                        if('*' == quality[i]) //skip read with no quality
                        {
                            i++; j++;
                            continue;
                        }
                        //TODO use Map instead
                        switch (read[i]){
                            case 'A':
                                totalQuality[0][j][0] += quality[i];
                                totalQuality[0][j][1]++;
                                break;
                            case 'C':
                                totalQuality[1][j][0] += quality[i];
                                totalQuality[1][j][1]++;
                                break;
                            case 'G':
                                totalQuality[2][j][0] += quality[i];
                                totalQuality[2][j][1]++;
                                break;
                            case 'T':
                                totalQuality[3][j][0] += quality[i];
                                totalQuality[3][j][1]++;
                                break;
                            case '-':
                                totalQuality[4][j][0] += quality[i];
                                totalQuality[4][j][1]++;
                                break;
                        }
                        i++; j++;
                    }
                }

                logger.info("Done. Start constructing the result quality string..");

                //print output quality string
                PrintWriter out = new PrintWriter(new FileWriter("contig_" + bedRecord.getStartIndex() + "-"
                + bedRecord.getStopIndex() + ".fastq"));

                //TODO expand with new symbol for DEL
                char[] dictNucleotid = {'A', 'C', 'G', 'T', 'D'};
                //restore the reference genome: find the most frequent symbol in each position. In most cases it is enough
                char[] reference = new char[bedRecord.getStopIndex() - bedRecord.getStartIndex()];
                for(int i = 0; i < reference.length; ++i)
                {
                    int max_count = totalQuality[0][i][1], nucleotid_index = 0;
                    for(int k = 1; k < dictNucleotid.length; ++k)
                        if(totalQuality[k][i][1] > max_count)
                        {
                            max_count = totalQuality[k][i][1];
                            nucleotid_index = k;
                        }

                    if(max_count > 0)
                        reference[i] = dictNucleotid[nucleotid_index];
                    else
                        reference[i] = 'N';
                }

                //print 4 different quality strings: first - is the quality if A is in every position, second - if C, and so on
                for(int k = 0; k < dictNucleotid.length; ++k){
                    out.println("@BED:" + bedRecord.getStartIndex() + " - " + bedRecord.getStopIndex()+ "(" + dictNucleotid [k] + ")");
                    out.println(reference);
                    out.println("+");
                    for(int j = 0; j < bedRecord.getStopIndex() - bedRecord.getStartIndex(); ++j)
                        if(totalQuality[k][j][1] > 0)
                            out.print((char)(QualityExtractor.QUAL_BASE + totalQuality[k][j][0] / totalQuality[k][j][1]));
                        else
                            out.print('*'); //no quality is stored
                    out.print("\n");
                }
                out.close();
            }
        }catch (Exception e){
            System.out.println("Error while data processing: " + e.toString());
        }
    }
    public static final byte QUAL_BASE = 33;
}
