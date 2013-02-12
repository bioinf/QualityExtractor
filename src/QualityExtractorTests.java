import difflib.DiffUtils;
import difflib.Patch;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 03.01.13
 * Time: 13:51
 * To change this template use File | Settings | File Templates.
 */
public class QualityExtractorTests {

    private static List<String> fileToLines(String filename) {
        List<String> lines = new LinkedList<String>();
        String line;
        try {
            BufferedReader in = new BufferedReader(new FileReader(filename));
            while ((line = in.readLine()) != null) {
                lines.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return lines;
    }

    private static String inputBed = "./TestData/cftr.bed";
    private static String genomeRef = "/Users/Kos/Dropbox/Bioinf/Data/chromosomes/hg19.fa";

    @Before
    public void setUp() throws Exception {

    }

    @After
    public void tearDown() throws Exception {

    }

    @Test
    public void testLeftOverlapFASTQ() throws Exception {
        String inputBam = "./TestData/test_left_overlap.sam.sorted.bam";
        QualityExtractor.createBamIndex(new File(inputBam), new File(inputBam + ".bai"));
        QualityExtractor.calcMeanQuality(QualityExtractor.extractReadsFromBAM(new File(inputBam), new File(inputBed)),
                genomeRef, QualityExtractor.OutputFormat.FASTQ, ".");

        List<String> original = fileToLines("./TestData/test_left_overlap_117138345-117138546.fastq");
        List<String> testResult = fileToLines("./region_117138345-117138546.fastq");
        Patch patch = DiffUtils.diff(original, testResult);
        assertEquals(patch.getDeltas().size(), 0);
    }

    @Test
    public void testLeftOverlapVCF() throws Exception {
        String inputBam = "./TestData/test_left_overlap.sam.sorted.bam";
        QualityExtractor.createBamIndex(new File(inputBam), new File(inputBam + ".bai"));
        QualityExtractor.calcMeanQuality(QualityExtractor.extractReadsFromBAM(new File(inputBam), new File(inputBed)),
                genomeRef, QualityExtractor.OutputFormat.VCF, ".");
        List<String> original = fileToLines("./TestData/test_left_overlap_117138345-117138546.vcf");
        List<String> testResult = fileToLines("./region_117138345-117138546.vcf");
        Patch patch = DiffUtils.diff(original, testResult);
        assertEquals(patch.getDeltas().size(), 0);
    }

    @Test
    public void testRightOverlapFASTQ() throws Exception {
        String inputBam = "./TestData/test_right_overlap.sam.sorted.bam";
        QualityExtractor.calcMeanQuality(QualityExtractor.extractReadsFromBAM(new File(inputBam), new File(inputBed)),
                genomeRef, QualityExtractor.OutputFormat.FASTQ, ".");
        List<String> original = fileToLines("./TestData/test_right_overlap_117138345-117138546.fastq");
        List<String> testResult = fileToLines("./region_117138345-117138546.fastq");
        Patch patch = DiffUtils.diff(original, testResult);
        assertEquals(patch.getDeltas().size(), 0);
    }

    @Test
    public void testRightOverlapVCF() throws Exception {
        String inputBam = "./TestData/test_right_overlap.sam.sorted.bam";
        QualityExtractor.calcMeanQuality(QualityExtractor.extractReadsFromBAM(new File(inputBam), new File(inputBed)),
                genomeRef, QualityExtractor.OutputFormat.VCF, ".");
        List<String> original = fileToLines("./TestData/test_right_overlap_117138345-117138546.vcf");
        List<String> testResult = fileToLines("./region_117138345-117138546.vcf");
        Patch patch = DiffUtils.diff(original, testResult);
        assertEquals(patch.getDeltas().size(), 0);
    }

    @Test
    public void testBothOverlapFASTQ() throws Exception {
        String inputBam = "./TestData/test_both_overlap.sam.sorted.bam";
        QualityExtractor.calcMeanQuality(QualityExtractor.extractReadsFromBAM(new File(inputBam), new File(inputBed)),
                genomeRef, QualityExtractor.OutputFormat.FASTQ, ".");
        List<String> original = fileToLines("./TestData/test_both_overlap_117138345-117138546.fastq");
        List<String> testResult = fileToLines("./region_117138345-117138546.fastq");
        Patch patch = DiffUtils.diff(original, testResult);
        assertEquals(patch.getDeltas().size(), 0);
    }

    @Test
    public void testBothOverlapVCF() throws Exception {
        String inputBam = "./TestData/test_both_overlap.sam.sorted.bam";
        QualityExtractor.calcMeanQuality(QualityExtractor.extractReadsFromBAM(new File(inputBam), new File(inputBed)),
                genomeRef, QualityExtractor.OutputFormat.VCF, ".");
        List<String> original = fileToLines("./TestData/test_both_overlap_117138345-117138546.vcf");
        List<String> testResult = fileToLines("./region_117138345-117138546.vcf");
        Patch patch = DiffUtils.diff(original, testResult);
        assertEquals(patch.getDeltas().size(), 0);
    }

}
