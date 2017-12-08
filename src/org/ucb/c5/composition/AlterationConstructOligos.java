package org.ucb.c5.composition;

import javafx.util.Pair;
import org.ucb.c5.utils.RevComp;

import java.util.ArrayList;

/**
 * The protocol for constructing promoter alteration constructs for each TF-
 * motif location in the promoter and the UTR of the specified gene relies on
 * an overlap-extension PCR approach that removes an entire motif.
 *
 * @author Manraj Gill
 */
public class AlterationConstructOligos {

    private RevComp rc;

    public void initiate() throws Exception {
        rc = new RevComp();
        rc.initiate();
    }

    public ArrayList<Pair<String, ArrayList<String>>> run(ArrayList<Pair<String, ArrayList<Integer>>> TFmotifs, String promoterAndFivePrimeUTR) throws Exception {
        ArrayList<Pair<String, ArrayList<String>>> alterationConstructOligos = new ArrayList<>();

        for (Pair<String, ArrayList<Integer>> TFmotif : TFmotifs) {
            String Motif_ID = TFmotif.getKey();
            ArrayList<Integer> motifLocation = TFmotif.getValue();
            int motifLength = motifLocation.get(0);
            for (int i = 1; i < motifLocation.size(); i++) {
                int locationStartIndex = motifLocation.get(i);
                // Can only construct primers for sites that are at least 15 basepairs from the edges
                if (locationStartIndex < 15 | locationStartIndex > (promoterAndFivePrimeUTR.length() - motifLength - 16)) {
                    continue;
                }

                // Construct the oligo with left and right components that exclude the motif
                int startingIndexOfOligoLeft = locationStartIndex - 15;
                int endingIndexOfOligoLeft = locationStartIndex;
                int startingIndexOfOligoRight = locationStartIndex + motifLength;
                int endingIndexOfOligoRight = startingIndexOfOligoRight + 15;
                String oligoLeft = promoterAndFivePrimeUTR.substring(startingIndexOfOligoLeft, endingIndexOfOligoLeft);
                String oligoRight = promoterAndFivePrimeUTR.substring(startingIndexOfOligoRight, endingIndexOfOligoRight);
                String forwardOligo = oligoLeft.concat(oligoRight);

                // Derive the reverse complement of the oligo that will serve as a reverse primer
                String reverseOligo = rc.run(forwardOligo);

                // Add to the arraylist named 'alterationConstructOligos' the forward and reverse oligos
                // for a key corresponding to Motif_ID and the site on the sequence
                String appendedMotif_ID = Motif_ID.concat("_");
                appendedMotif_ID = appendedMotif_ID.concat(Integer.toString(locationStartIndex));
                ArrayList<String> oligos = new ArrayList<>();
                oligos.add(0, forwardOligo);
                oligos.add(1, reverseOligo);
                Pair<String, ArrayList<String>> alterationOligos = new Pair<>(appendedMotif_ID, oligos);
                alterationConstructOligos.add(alterationOligos);
            }
        }

        return alterationConstructOligos;
    }
}
