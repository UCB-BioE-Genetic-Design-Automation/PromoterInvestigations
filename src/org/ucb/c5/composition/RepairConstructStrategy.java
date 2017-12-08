package org.ucb.c5.composition;

import org.ucb.c5.utils.RevComp;
import java.util.ArrayList;
import java.util.HashMap;

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
 * Initiating RepairConstructStrategy constructs a library of the restriction
 * endonuclease sites present in the Multiple Cloning Site (MCS) of the pCR2.1-
 * TOPO plasmid that are sticky-end and unique 6+ or 8+ cutters.
 *
 * Running RepairConstructStrategy determines which restriction endonucleases
 * can be used for cloning the entire locus into the vector backbone. And it
 * then constructs primers with using these sites.
 *
 * @author Manraj Gill
 */
public class RepairConstructStrategy {

    private RevComp rc;
    private ArrayList<String> RestrictionEndonucleaseSites;
    private HashMap<String, String> RestrictionEndonucleases;

    public void initiate() throws Exception {
        rc = new RevComp();
        rc.initiate();

        RestrictionEndonucleaseSites = new ArrayList<>();
        RestrictionEndonucleases = new HashMap<>();

        RestrictionEndonucleaseSites.add("AAGCTT");
        RestrictionEndonucleases.put("AAGCTT", "HindIII");
        RestrictionEndonucleaseSites.add("GGATCC");
        RestrictionEndonucleases.put("GGATCC", "BamHI");
        RestrictionEndonucleaseSites.add("GGTACC");
        RestrictionEndonucleases.put("GGTACC", "KpnI");
        RestrictionEndonucleaseSites.add("GAGCTC");
        RestrictionEndonucleases.put("GAGCTC", "SacI");
        RestrictionEndonucleaseSites.add("ACTAGT");
        RestrictionEndonucleases.put("ACTAGC", "SpeI");
        RestrictionEndonucleaseSites.add("CTTAAG");
        RestrictionEndonucleases.put("CTTAAG", "AflII");
        RestrictionEndonucleaseSites.add("GCGGCCGC");
        RestrictionEndonucleases.put("GCGGCCGC", "NotI");
        RestrictionEndonucleaseSites.add("CTCGAG");
        RestrictionEndonucleases.put("CTCGAG", "XhoI");
        RestrictionEndonucleaseSites.add("CTCGAG");
        RestrictionEndonucleases.put("CTCGAG", "NsiI");
        RestrictionEndonucleaseSites.add("TCTAGA");
        RestrictionEndonucleases.put("TCTAGA", "XbaI");
        RestrictionEndonucleaseSites.add("GGGCCC");
        RestrictionEndonucleases.put("GGGCCC", "ApaI");
    }

    public ArrayList<String> run(String LeftHomologyArm, String PromoterAndFivePrimeUTR, String RightHomologyArm) throws Exception {
        // Construct a string that combines the three elements of the locus
        String entireLocus = LeftHomologyArm.concat(PromoterAndFivePrimeUTR);
        entireLocus = entireLocus.concat(RightHomologyArm);
        String entireLocusReverseComplement = rc.run(entireLocus);

        // Obtain the sites of the restriction endonucleases that are not present in the given gene locus
        String restrictionEndonucleaseSite1 = null;
        String restrictionEndonucleaseSite2 = null;
        int count = 0;
        for (String restrictionEndonucleaseSite : RestrictionEndonucleaseSites) {
            if (count == 2) {
                break;
            }
            if (!entireLocus.contains(restrictionEndonucleaseSite)) {
                if (count == 0) {
                    restrictionEndonucleaseSite1 = restrictionEndonucleaseSite;
                } else {
                    restrictionEndonucleaseSite2 = restrictionEndonucleaseSite;
                }
                count += 1;
            }
        }
        if (restrictionEndonucleaseSite1 == null || restrictionEndonucleaseSite2 == null) {
            throw new Exception("Could not devise a cloning strategy for this gene locus because pCR2.1-TOPO MCS does not contain unique sites for this locus!");
        }
        // Obtain the names of the restriction endonculeases corresponding to these sites that are not present in the given gene locus
        String restrictionEndonuclease1 = RestrictionEndonucleases.get(restrictionEndonucleaseSite1);
        String restrictionEndonuclease2 = RestrictionEndonucleases.get(restrictionEndonucleaseSite2);

        // Construct primers by appending initially a 5' buffer and 18bp of the edge sequences to the restriction endonuclease sites
        String forwardPrimer = "ATAT";
        forwardPrimer = forwardPrimer.concat(restrictionEndonucleaseSite1);
        forwardPrimer = forwardPrimer.concat(entireLocus.substring(0, 18));
        String reversePrimer = "ATAT";
        reversePrimer = reversePrimer.concat(restrictionEndonucleaseSite2);
        reversePrimer = reversePrimer.concat(entireLocusReverseComplement.substring(0, 18));

        // Add the elements to return to an arraylist named repairConstructStrategy
        ArrayList<String> repairConstructStrategy = new ArrayList<>();
        repairConstructStrategy.add(0, restrictionEndonuclease1);
        repairConstructStrategy.add(1, restrictionEndonuclease2);
        repairConstructStrategy.add(2, forwardPrimer);
        repairConstructStrategy.add(3, reversePrimer);
        repairConstructStrategy.add(4, Integer.toString(entireLocus.length()) + "bp");

        return repairConstructStrategy;
    }
}
