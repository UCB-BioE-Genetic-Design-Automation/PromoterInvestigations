package org.ucb.c5.composition;

import javafx.util.Pair;
import org.ucb.c5.utils.RevComp;
import java.util.ArrayList;

/**
 * The vector backbone into which the repair construct will be synthesized for
 * repair at the editing locus will be from the 'TOPO TA Cloning Kit with PCR2.1'
 * from Thermo Fisher Scientific
 *
 * Plasmid: pCR2.1-TOPO
 * Stored in this package as an .rtf file and a .dna file
 *
 * URL: https://www.addgene.org/vector-database/2285/
 *
 *
 * The protocol for constructing repair constructs for each transcription factor
 * motif location in the promoter and the UTR of the specified gene relies on
 * an overlap-extension PCR approach that removes an entire motif.
 *
 * @author Manraj Gill
 */
public class RepairConstructOligos {

    private RevComp rc;

    public void initiate() throws Exception {
        rc = new RevComp();
        rc.initiate();
    }

    public ArrayList<Pair<String, ArrayList<String>>> run(ArrayList<Pair<String, ArrayList<Integer>>> TFmotifs, String promoterAndFivePrimeUTR) throws Exception {
        ArrayList<Pair<String, ArrayList<String>>> repairOligos = new ArrayList<>();
        return repairOligos;
    }
}
