package org.ucb.c5.composition;

import java.util.ArrayList;
import java.util.HashMap;
import org.ucb.c5.utils.FileUtils;

/**
 * Data of genomic sequence from:
 *
 * "Promoter + 5' UTR.txt" file is an ENSEMBLE data extraction from the dataset of
 * Human genes (GRCh38.p10) filtered to only include genes with NCBI gene IDs.
 * The attributes in this file are its Gene stable ID, Gene name, Transcript
 * stable ID and its 5' UTR sequence with 1500bp upstream sequence.
 *
 * URL: http://www.ensembl.org/downloads.html
 *
 * Initiation of org.ucb.c5.composition.UpstreamGenomicLocus parses the data from the database
 * Running of org.ucb.c5.composition.UpstreamGenomicLocus outputs the LHA [0th index of returned
 * ArrayList and the promoter + 5' UTR sequence [1st index of the ArrayList] for the input gene.
 *
 * @author Manraj Gill
 */
public class UpstreamGenomicLocus {

    private HashMap<String, String> LHAs;
    private HashMap<String, String> promotersAndFivePrimeUTRs;

    public void initiate() throws Exception {
        String PromoterAndUTRinformation = FileUtils.readResourceFile("Promoter + 5' UTR.txt");
        String[] genes = PromoterAndUTRinformation.split(">");
        LHAs = new HashMap<>();
        promotersAndFivePrimeUTRs = new HashMap<>();
        for (int i = 1; i < genes.length; i++) {
            String gene = genes[i];
            String[] lines = gene.split("\\r|\\r?\\n");
            // For the cases: 'Sequence unavailable', skip without
            // including the information in the HashMaps.
            if (lines.length < 3) {
                continue;
            }
            String geneInformation = lines[0];
            // Removes the Gene stable ID from geneInformation
            String geneName = geneInformation.substring(16);
            // Removes the Transcript stable ID from geneInformation
            int endOfName = geneName.length() - 16;
            geneName = geneName.substring(0, endOfName);
            String complete = "";
            for (int j = 1; j < lines.length; j++) {
                complete += lines[j];
            }
            if (complete.length() < 1500) {
                continue;
            }
            if (LHAs.containsKey(geneName)) {
                continue;
            }
            String LHA = complete.substring(0, 500);
            LHAs.put(geneName, LHA);
            String PromoterAndFivePrimeUTR = complete.substring(500);
            promotersAndFivePrimeUTRs.put(geneName, PromoterAndFivePrimeUTR);
        }
    }

    public ArrayList<String> run(String geneName) throws Exception {
        ArrayList<String> toReturn = new ArrayList<>(2);
        if (LHAs.containsKey(geneName) && promotersAndFivePrimeUTRs.containsKey(geneName)) {
            toReturn.add(0, LHAs.get(geneName));
            toReturn.add(1, promotersAndFivePrimeUTRs.get(geneName));
        }
        else {
            throw new IllegalArgumentException("The input gene name is not valid.");
        }
        return toReturn;
    }

    public static void main(String[] args) throws Exception {
        UpstreamGenomicLocus ugl = new UpstreamGenomicLocus();
        ugl.initiate();
        ArrayList<String> TERTugl = ugl.run("TERT");
        System.out.println(TERTugl.get(1));
    }
}
