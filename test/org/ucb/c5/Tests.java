package org.ucb.c5;

import org.junit.Test;
import org.junit.BeforeClass;
import static org.junit.Assert.*;
import org.ucb.c5.composition.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import javafx.util.Pair;

/**
 * Tests PromoterInvestigations functionality by comparing outputs
 * to human generated oligos and strategies for the hTERT locus.
 *
 * hTERT is the human Telomerase Reverse Transcriptase gene and this
 * gene is transcriptionally regulated during both cellular differen-
 * tiation and tumorigenesis by being silenced and re-activated, res-
 * pectively. Therefore, its promoter has been the subject of many
 * investigations aimed at understanding the mechanism of action of
 * the transcription factors that bind.
 *
 * In these tests, the known binding site of a transcription factor
 * Myc at the TERT locus ~ 190 base pairs upstream of the Trancript-
 * ional Start Site is confirmed. Myc binds enhancer box sequences
 * (E-boxes) and recruits histone acetyltranferases to activate the
 * trancription of genes.
 *
 * @author Manraj Gill
 */
public class Tests {

    public Tests() {}

    private static DownstreamGenomicLocus dgl;
    private static UpstreamGenomicLocus ugl;
    private static TFmotifs tfms;
    private static Cas9ConstructOligos c9co;
    private static RepairConstructStrategy rcs;

    @BeforeClass
    public static void testClassSetUp() throws Exception {
        dgl = new DownstreamGenomicLocus();
        ugl = new UpstreamGenomicLocus();
        tfms = new TFmotifs();
        c9co = new Cas9ConstructOligos();
        rcs = new RepairConstructStrategy();
        dgl.initiate();
        ugl.initiate();
        tfms.initiate();
        c9co.initiate();
        rcs.initiate();
    }

    @Test
    public void DownstreamGenomicLocusTest() throws Exception {
        String RightHomologyArm = dgl.run("TERT");
        String first500bpAfterATGofTERT = "atgccgcgcgctccccgctgccgagccgtgcgctccctgctgcgcagccactaccgcgaggtgctgccgctggccacgttcgtgcggcgcctggggccccagggctggcggctggtgcagcgcggggacccggcggctttccgcgcgctggtggcccagtgcctggtgtgcgtgccctgggacgcacggccgccccccgccgccccctccttccgccaggtgggcctccccggggtcggcgtccggctggggttgagggcggccggggggaaccagcgacatgcggagagcagcgcaggcgactcagggcgcttcccccgcaggtgtcctgcctgaaggagctggtggcccgagtgctgcagaggctgtgcgagcgcggcgcgaagaacgtgctggccttcggcttcgcgctgctggacggggcccgcgggggcccccccgaggccttcaccaccagcgtgcgcagctacctgcccaacacggtgaccgacgcactgcgg";
        assertEquals(RightHomologyArm, first500bpAfterATGofTERT.toUpperCase());
    }

    @Test
    public void UpstreamGenomicLocusTest() throws Exception {
        ArrayList<String> upstreamGenomicLocus = ugl.run("TERT");
        String promoterAnd5primeUTRofTERT = "agtggccgtgtggcttctactgctgggctggaagtcgggcctcctagctctgcagtccgaggcttggagccaggtgcctggaccccgaggttgccctccaccctgtgcgggcgggatgtgaccagatgttggcctcatctgccagacagagtgccggggcccagggtcaaggccgttgtggctggtgtgaggcgcccggtgcgcggccagcaggagcgcctggctccatttcccaccctttctcgacgggaccgccccggtgggtgattaacagatttggggtggtttgctcatggtggggacccctcgccgcctgagaacctgcaaagagaaatgacgggcctgtgtcaaggagcccaagtcgcggggaagtgttgcagggaggcactccgggaggtcccgcgtgcccgtccagggagcaatgcgtcctcgggttcgtccccagccgcgtctacgcgcctccgtcctccccttcacgtccggcattcgtggtgcccggagcccgacgccccgcgtccggacctggaggcagccctgggtctccggatcaggccagcggccaaagggtcgccgcacgcacctgttcccagggcctccacatcatggcccctccctcgggttaccccacagcctaggccgattcgacctctctccgctggggccctcgctggcgtccctgcaccctgggagcgcgagcggcgcgcgggcggggaagcgcggcccagacccccgggtccgcccggagcagctgcgctgtcggggccaggccgggctcccagtggattcgcgggcacagacgcccaggaccgcgcttcccacgtggcggagggactggggacccgggcacccgtcctgccccttcaccttccagctccgcctcctccgcgcggaccccgccccgtcccgacccctcccgggtccccggcccagccccctccgggccctcccagcccctccccttcctttccgcggccccgccctctcctcgcggcgcgagtttcaggcagcgctgcgtcctgctgcgcacgtgggaagccctggccccggccacccccgcg";
        String upstreamOfPromoter500bp = "taaaattgtgttttctatgttggcttctctgcagagaaccagtgtaagctacaacttaacttttgttggaacaaattttccaaaccgcccctttgccctagtggcagagacaattcacaaacacagccctttaaaaaggcttagggatcactaaggggatttctagaagagcgacctgtaatcctaagtatttacaagacgaggctaacctccagcgagcgtgacagcccagggagggtgcgaggcctgttcaaatgctagctccataaataaagcaatttcctccggcagtttctgaaagtaggaaaggttacatttaaggttgcgtttgttagcatttcagtgtttgccgacctcagctacagcatccctgcaaggcctcgggagacccagaagtttctcgccccttagatccaaacttgagcaacccggagtctggattcctgggaagtcctcagctgtcctgcggttgtgccggggccccaggtctggaggggacc";
        String LeftHomologyArm = upstreamGenomicLocus.get(0);
        assertEquals(LeftHomologyArm, upstreamOfPromoter500bp.toUpperCase());
        String uglPromoterAnd5primeUTR = upstreamGenomicLocus.get(1);
        assertEquals(uglPromoterAnd5primeUTR, promoterAnd5primeUTRofTERT.toUpperCase());
    }

    @Test
    public void TFmotifsTest() throws Exception {
        ArrayList<String> upstreamGenomicLocus = ugl.run("TERT");
        String PromoterAndFivePrimeUTR = upstreamGenomicLocus.get(1);
        ArrayList<ArrayList> TFmotifsInSequence = tfms.run(PromoterAndFivePrimeUTR);
        ArrayList<HashMap> TFinformation = TFmotifsInSequence.get(1);
        HashMap<String, String> tfNames = TFinformation.get(0);
        boolean MycPresent = false;
        for (Map.Entry<String, String> entry : tfNames.entrySet()) {
            String tfName = entry.getValue();
            if (tfName.equals("MYC")) {
                MycPresent = true;
                break;
            }
        }
        assertTrue(MycPresent);
    }

    @Test
    public void Cas9ConstructOligosTest() throws Exception {
        ArrayList<String> upstreamGenomicLocus = ugl.run("TERT");
        String PromoterAndFivePrimeUTR = upstreamGenomicLocus.get(1);
        ArrayList<ArrayList> TFmotifsInSequence = tfms.run(PromoterAndFivePrimeUTR);
        ArrayList<Pair<String, ArrayList<Integer>>> TranscriptionFactors = TFmotifsInSequence.get(0);
        ArrayList<Pair<String, ArrayList<String>>> sgRNAs = c9co.run(TranscriptionFactors, PromoterAndFivePrimeUTR);
        HashMap<String, ArrayList<String>> hashMapOfsgRNAs = new HashMap<>();
        for (Pair<String, ArrayList<String>> sgRNA : sgRNAs) {
            String appendedMotif_ID = sgRNA.getKey();
            ArrayList<String> sgRNAoligos = sgRNA.getValue();
            hashMapOfsgRNAs.put(appendedMotif_ID, sgRNAoligos);
        }
        // The Myc motif recognized at -186 base pairs from the TSS has a motif_id of M4610_1.02
        // -186bp from the TSS corresponds to 814 basepairs from the start of the promotoer (defined as a 1000bp segment)
        ArrayList<String> TERTsgRNAoligos = hashMapOfsgRNAs.get("M4610_1.02_814");
        String sgUP = TERTsgRNAoligos.get(0);
        String sgDown = TERTsgRNAoligos.get(1);
        String expectedSgUP = "caccgCCAGGACCGCGCTTCCCACG";
        String expectedSgDown = "aaacCGTGGGAAGCGCGGTCCTGGc";
        assertEquals(sgUP, expectedSgUP.toUpperCase());
        assertEquals(sgDown, expectedSgDown.toUpperCase());
    }

    @Test
    public void AlterationConstructOligosTest() throws Exception {
        ArrayList<String> upstream = ugl.run("TERT");
        String RightHomologyArm = dgl.run("TERT");
        String LeftHomologyArm = upstream.get(0);
        String PromoterAndFivePrimeUTR = upstream.get(1);
        ArrayList<String> repairStrategy = rcs.run(LeftHomologyArm, PromoterAndFivePrimeUTR, RightHomologyArm);
        // The TERT locus does not contain a HindIII or a BamHI site and given that the algorithm is designed
        // to pick the first restriction endonucleases that aren't present in the locus, we deterministically
        // obtain a repair construct strategy that employs HindIII and BamHI for this gene.
        assertEquals(repairStrategy.get(0), "HindIII");
        assertEquals(repairStrategy.get(1), "BamHI");
        // The expected size of the amplicon should equal the combined length of the LHA, promoter, 5' UTR and RHA
        int expectedSizeOfAmplicon = Integer.parseInt(repairStrategy.get(4).substring(0, 4));
        assertEquals(expectedSizeOfAmplicon, RightHomologyArm.length() + LeftHomologyArm.length() + PromoterAndFivePrimeUTR.length());

        String expectedForward = "atataagctttaaaattgtgttttctat";
        String expectedReverse = "atatggatccccgcagtgcgtcggtcac";
        assertEquals(repairStrategy.get(2), expectedForward.toUpperCase());
        assertEquals(repairStrategy.get(3), expectedReverse.toUpperCase());
    }
}
