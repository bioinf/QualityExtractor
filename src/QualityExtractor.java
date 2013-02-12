import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
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

        logger.info("Launch Arguments count: " + args.length + "\nArguments list: ");
        for(String str : args){
            logger.info(str + ";");
        }
        if(args.length != 4) {
            logger.error("Not enough arguments to launch. Quit.");
            return;
        }

        logger.info("Start processing...");
        createBamIndex(new File(args[0]), new File(args[0] + ".bai"));
        calcMeanQuality(extractReadsFromBAM(new File(args[0]), new File(args[1])),
                args[2], OutputFormat.VCF, args[3]);
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
            if(outputBamIndexFile.exists()){
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
            final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile, new File(inputSamOrBamFile.getAbsolutePath() + ".bai"));
            if(!inputSam.hasIndex())
                throw new Exception("Index file is invalid");
            while (inputBed.hasNext()){
                BEDRecord record = inputBed.getBEDRecord();
                logger.info("BEDRecord: " + record.getContigName()+ "\t" + record.getStartIndex() + '\t' + record.getStopIndex());
                SAMRecordIterator it = inputSam.query(record.getContigName(), record.getStartIndex(), record.getStopIndex(), false);
                List<SAMRecord> list = new ArrayList<SAMRecord>();
                int totalRecords = 0;
                while (it.hasNext()){
                    SAMRecord rec = it.next();
                    list.add(rec);
                    totalRecords++;
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
    public static void getAlignedQualities(final SAMRecord record, byte[] outputQuality, byte[] outputRead){
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
    public static void calcMeanQuality(final Map<BEDRecord, List<SAMRecord>> dict, String refGenome, OutputFormat format, String folder) {
        Logger logger = Logger.getLogger(QualityExtractor.class);
        IndexedFastaSequenceFile hg = null;
        PrintWriter out = null;

        try{
            Set s = dict.entrySet();
            Iterator it = s.iterator();
            //TODO obtain path from environment variables
            hg = new IndexedFastaSequenceFile(new File(refGenome));

            //Iterate through map records
            while(it.hasNext()){
                Map.Entry m =(Map.Entry)it.next();
                BEDRecord bedRecord = (BEDRecord) m.getKey();
                List<SAMRecord> records = (List<SAMRecord>) m.getValue();
                logger.info("Parsing reference sequence: " + bedRecord.getContigName() + " " + bedRecord.getStartIndex() +
                        " " + bedRecord.getStopIndex() + " with " + records.size() + " sequences");

                Map<String, ArrayList<Pair>> dataSet = new LinkedHashMap<String, ArrayList<Pair>>();
                for(char nucleotid : dictNucleotid){
                    ArrayList<Pair> list = new ArrayList<Pair>();
                    for(int i = 0; i < bedRecord.getStopIndex() - bedRecord.getStartIndex() + 1; ++i){
                        list.add(i, new Pair(0, 0));
                    }
                    dataSet.put(String.valueOf(nucleotid), list);
                }

                logger.info("Start calculating the mean value for quality...");
                for(SAMRecord samRecord : records){
                    if(0 == samRecord.getAlignmentStart()) {
                        continue;     //unaligned read. Skip
                    }

                    byte[] quality = new byte[samRecord.getAlignmentEnd() - samRecord.getAlignmentStart() + 1];
                    byte[] read = new byte[samRecord.getAlignmentEnd() - samRecord.getAlignmentStart() + 1];
                    QualityExtractor.getAlignedQualities(samRecord, quality, read);

                    AlignmentIndicies indicies = getIndicies(samRecord, bedRecord);

                    //find indicies to correctly copy overlapped record
                    int i = indicies.samStart, // start index in SAMRecord
                            j = indicies.bedStart, // start index in BEDRecord
                            copySize = indicies.copySize; //num of bytes to copy from SAMRecord

                    //copy SAMRecord to BEDRecord boundaries
                    copySize += i; //if i != 0, because we need to iterate copySize times
                    while(i < copySize){
                        if('*' == quality[i]){ //skip read with no quality
                            i++; j++;
                            continue;
                        }
                        Pair currentPair = dataSet.get(String.valueOf((char)read[i])).get(j);
                        currentPair.first += quality[i];
                        currentPair.second += 1;
                        dataSet.get(String.valueOf((char)read[i])).set(j, currentPair);
                        i++; j++;
                    }
                }

                logger.info("Done. Start constructing the result quality string..");

                //print output quality string
                out = new PrintWriter(new FileWriter(folder + "/region_" + bedRecord.getStartIndex() + "-"
                + bedRecord.getStopIndex() + (format.equals(OutputFormat.FASTQ) ? ".fastq" : ".vcf")));

                ReferenceSequence refSequence =  hg.getSubsequenceAt(bedRecord.getContigName(), bedRecord.getStartIndex(), bedRecord.getStopIndex());
                switch (format){
                    case VCF:
                        printVcf(out, dataSet, bedRecord, refSequence);
                        break;
                    case FASTQ:
                        printFastq(out, dataSet, bedRecord, refSequence);
                        break;
                }
                out.close();
            }
        } catch (Exception e){
            if(out != null)
                out.close();
            System.out.println("Error while data processing: " + e.toString());
        }
    }

    private static void printFastq(PrintWriter out, final Map<String, ArrayList<Pair>> dataSet, final BEDRecord bedRecord, final ReferenceSequence refSequence){
        //print 4 different quality strings: first - is the quality if A is in every position, second - if C, and so on
        for(char nucleotid : dictNucleotid){
            out.println("@BED:" + bedRecord.getStartIndex() + " - " + bedRecord.getStopIndex()+ "(" + nucleotid + ")");
            out.println(new String(refSequence.getBases()));
            out.println("+");
            for(int j = 0; j < bedRecord.getStopIndex() - bedRecord.getStartIndex(); ++j){
                if(dataSet.get(String.valueOf(nucleotid)).get(j).second > 0){
                    out.print((char)(QualityExtractor.QUAL_BASE + dataSet.get(String.valueOf(nucleotid)).get(j).first / dataSet.get(String.valueOf(nucleotid)).get(j).second));
                } else {
                    out.print('*'); //no quality is stored
                }
            }
            out.print("\n");
        }
    }

    private static void printVcf(PrintWriter out, final Map<String, ArrayList<Pair>> dataSet, final BEDRecord bedRecord, final ReferenceSequence refSequence){
        out.println("##fileformat=VCFv4.0");
        out.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        byte[] ref = refSequence.getBases();
        for(int i = 0; i < bedRecord.getStopIndex() - bedRecord.getStartIndex() + 1; ++i){
            char bestNucleotid = 'N';
            int bestNucleotidCount = 0;
            int bestNucleotidQuality = 0;
            for(char nucleotid : dictNucleotid){
                Pair currentPair = dataSet.get(String.valueOf(nucleotid)).get(i);
                if(currentPair.second > bestNucleotidCount){
                    bestNucleotid = nucleotid;
                    bestNucleotidQuality =  currentPair.first;
                    bestNucleotidCount = currentPair.second;
                }
            }
            StringBuilder str = new StringBuilder();
            str.append(bedRecord.getContigName() + "\t" +       //CHROM
                    (i + bedRecord.getStartIndex()) + "\t" +    //POS
                    ".\t" +                                     //ID
                    (char)ref[i] + "\t" +                       //REF
                    (bestNucleotid == (char)ref[i] ? "." : String.valueOf(bestNucleotid)) + "\t" +  //ALT
                    ((bestNucleotidCount == 0) ? "." : (bestNucleotidQuality / bestNucleotidCount)) + "\t" + //QUAL
                    ".\t" +                                     //FILTER
                    ".\t"                                       //INFO
            );
            out.println(str.toString());
        }
    }

    private static class Pair{

        Pair(int first, int second){
            this.first = first;
            this.second = second;
        }

        public int first;
        public int second;
    }

    /**
     * Calculate correct start indicies to fit record into BED region (SAM record can overlap BED region)
     * @param samRecord
     * @param bedRecord
     * @return start position in SAM sequence, start position in BED region and copy size
     */
    private final static AlignmentIndicies getIndicies(SAMRecord samRecord, BEDRecord bedRecord){
        int samStartPos = samRecord.getAlignmentStart(), samStopPos = samRecord.getAlignmentEnd(),
                bedStartPos = bedRecord.getStartIndex(), bedStopPos = bedRecord.getStopIndex();

        //SAMRecord can overlap BEDRecord, so we need to be careful
        if(samStartPos >= bedStartPos && samStopPos <= bedStopPos){
            //SAMRecord is included in BED boundaries
            return new AlignmentIndicies(0, samStartPos - bedStartPos, samStopPos - samStartPos + 1);
        } else if(samStartPos < bedStartPos && samStopPos <= bedStopPos){
            //SAMRecord overlaps the left boundary
            return new AlignmentIndicies(bedStartPos - samStartPos, 0, samStopPos - bedStartPos + 1);
        } else if(samStartPos >= bedStartPos && samStopPos > bedStopPos){
            //SAMRecord overlaps the right boundary
            return new AlignmentIndicies(0, samStartPos - bedStartPos, bedStopPos - samStartPos + 1);
        } else if(samStartPos < bedStartPos && samStopPos > bedStopPos){
            //SAMRecord overlaps both boundaries
            return new AlignmentIndicies(bedStartPos - samStartPos, 0, bedStopPos - bedStartPos + 1);
        }
        return new AlignmentIndicies(0, 0, 0);
    }

    private final static class AlignmentIndicies{
        public final int copySize;
        public final int bedStart;
        public final int samStart;

         AlignmentIndicies(int samStart, int bedStart, int copySize) {
            this.copySize = copySize;
            this.bedStart = bedStart;
            this.samStart = samStart;
        }
    }

    private static char[] dictNucleotid = {'A', 'C', 'G', 'T', '-'};

    private static final byte QUAL_BASE = 33;

    public static enum OutputFormat{
        VCF,
        FASTQ
    }
}
