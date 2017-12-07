package org.ucb.c5.composition;

import javafx.util.Pair;
import org.ucb.c5.utils.RevComp;
import java.util.ArrayList;

/**
 * Design of CRISPR ( C lustered R egularly I nterspaced S hort P alindromic R epeats)
 * Cas9 single guide RNA oligos adapted from Zhang Lab CRISPR Protocol for implementat
 * -ion of this CRSIPR/Cas9 system in mammalian cells by co-expressing the S. pyogenes
 * Cas9 (SpCas9) nuclease along with the sgRNA.
 *
 * Plasmid: pX330-U6-Chimeric_BB-CBh-hSpCas9
 * Stored in this package as a .dna file visualizable using SnapGene Viewer software
 *
 * Reference:
 * Multiplex Genome Engineering using CRISPR/Cas Systems.
 * Cong L, Ran FA, Cox D, Lin S, Barretto R, Habib N, Hsu PD,
 * Wu X, Jiang W, Marraffini LA, Zhang F.
 * Science . 2013 Jan 3. DOI: 10.1126/science.1231143.
 *
 * URL: https://www.addgene.org/crispr/zhang/#spcas9
 *
 * @author Manraj Gill
 */
public class Cas9ConstructOligos {

    RevComp rc;

    public void initiate() throws Exception {
        rc = new RevComp();
        rc.initiate();
    }

    public ArrayList<Pair<String, ArrayList<String>>> run(ArrayList<Pair<String, ArrayList<Integer>>> TFmotifs, String promoterAndFivePrimeUTR) throws Exception {
        ArrayList<Pair<String, ArrayList<String>>> sgRNAs = new ArrayList<>();

        for (Pair<String, ArrayList<Integer>> TFmotif : TFmotifs) {
            String Motif_ID = TFmotif.getKey();
            ArrayList<Integer> motifLocation = TFmotif.getValue();
            for (int i = 1; i < motifLocation.size(); i++) {
                int locationStartIndex = motifLocation.get(i);
                if (locationStartIndex <= 22) {
                    continue;
                }
                int indexOfFirstDownstreamPAM = promoterAndFivePrimeUTR.indexOf("GG", locationStartIndex);
                if (indexOfFirstDownstreamPAM == -1) {
                    continue;
                }
                int startingIndexOfsgRNA = indexOfFirstDownstreamPAM - 21;
                int endingIndexOfsgRNA = indexOfFirstDownstreamPAM - 2;
                String sgRNA = promoterAndFivePrimeUTR.substring(startingIndexOfsgRNA, endingIndexOfsgRNA + 1);

                // Construct oligoUp by appending BbsI cut site
                String oligoUp = "";
                oligoUp = oligoUp.concat("CACCG");
                oligoUp = oligoUp.concat(sgRNA);

                // Construct oligoUp by appending BbsI cut site to the reverse complement of the sgRNA
                String sgRNArc = rc.run(sgRNA);
                String oligoDown = "";
                oligoDown = oligoDown.concat("AAAC");
                oligoDown = oligoDown.concat(sgRNArc);
                oligoDown = oligoDown.concat("C");

                // Add to the arraylist named 'sgRNAs' oligos up and down as values
                // for a key corresponding to Motif_ID and the site on the sequence
                String appendedMotif_ID = Motif_ID.concat("_");
                appendedMotif_ID = appendedMotif_ID.concat(Integer.toString(locationStartIndex));
                ArrayList<String> oligos = new ArrayList<>();
                oligos.add(0, oligoUp);
                oligos.add(1, oligoDown);
                Pair<String, ArrayList<String>> MotifOligos = new Pair<>(appendedMotif_ID, oligos);
                sgRNAs.add(MotifOligos);
            }
        }

        return sgRNAs;
    }

    public static void main(String[] args) throws Exception {
    }
}
