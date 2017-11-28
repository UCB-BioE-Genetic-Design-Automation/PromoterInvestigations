package org.ucb.c5.composition;

import javafx.util.Pair;

import java.util.ArrayList;
import java.util.HashMap;
import org.ucb.c5.utils.FileUtils;

/**
 * Data of genomic sequence from:
 *
 * "Exons + Introns" file is an ENSEMBLE data extraction from the dataset of Human
 * genes (GRCh38.p10) filtered to only include genes with NCBI gene IDs. The
 * attributes in this file are its Gene stable ID, Gene name, and 1000bp of its
 * unspliced genomic sequence downstream of the 5' UTR (also included).
 *
 * URL: http://www.ensembl.org/downloads.html
 *
 * Initiation of org.ucb.c5.composition.DownstreamGenomicLocus parses the data from the database
 * Running of org.ucb.c5.composition.DownstreamGenomicLocus outputs the RHA as a string for the
 * input gene.
 *
 * @author Manraj Gill
 */
public class DownstreamGenomicLocus {

    private HashMap<String, String> RHAs;

    public void initiate() throws Exception {
        String ExonAndIntronInformation = FileUtils.readResourceFile("Exons + Introns.txt");
        String[] allGenes = ExonAndIntronInformation.split(">");
        RHAs = new HashMap<>();
        for (int i = 1; i < allGenes.length; i++) {
            String gene = allGenes[i];
            String[] lines = gene.split("\\r|\\r?\\n");
            if (lines.length < 3) {
                continue;
            }
            String geneInformation = lines[0];
            String geneName = geneInformation.substring(16);
            String complete = "";
            for (int j = 1; j < lines.length; j++) {
                complete += lines[j];
            }
            // Need at least 500 basepairs of homology in the right homology arm (RHA)
            if (complete.length() < 500) {
                continue;
            }
            if (RHAs.containsKey(geneName)) {
                continue;
            }
            int lengthOfFivePrimeUTR = complete.length() - 1000;
            String RHA = complete.substring(lengthOfFivePrimeUTR, lengthOfFivePrimeUTR + 500);
            RHAs.put(geneName, RHA);
        }
    }

    public String run(String geneName) throws Exception {
        String toReturn = null;
        if (RHAs.containsKey(geneName)) {
            toReturn = RHAs.get(geneName);
        }
        else {
            throw new IllegalArgumentException("The input gene name is not valid.");
        }
        return toReturn;
    }

    public static void main(String[] args) throws Exception {
        DownstreamGenomicLocus dgl = new DownstreamGenomicLocus();
        dgl.initiate();
    }
}
