package org.ucb.c5.composition;

import javafx.util.Pair;

import java.util.ArrayList;

/**
 * @author Manraj Gill
 */
public class InvestigatePromoter {

    private DownstreamGenomicLocus dgl;
    private UpstreamGenomicLocus ugl;
    private TFmotifs tfms;
    private Cas9ConstructOligos c9co;
    private RepairConstructOligos rco;

    public void initiate() throws Exception {
        dgl = new DownstreamGenomicLocus();
        ugl = new UpstreamGenomicLocus();
        tfms = new TFmotifs();
        c9co = new Cas9ConstructOligos();
        rco = new RepairConstructOligos();
        dgl.initiate();
        ugl.initiate();
        tfms.initiate();
        c9co.initiate();
        rco.initiate();
    }

    public void run(String geneName) throws Exception {
        ArrayList<String> upstream = ugl.run(geneName);
        String RightHomologyArm = dgl.run(geneName);
        String LeftHomologyArm = upstream.get(0);
        String PromoterAndFivePrimeUTR = upstream.get(1);
        ArrayList<Pair<String, ArrayList<Integer>>> TranscriptionFactors = tfms.run(PromoterAndFivePrimeUTR);
        ArrayList<Pair<String, ArrayList<String>>> sgRNAs = c9co.run(TranscriptionFactors, PromoterAndFivePrimeUTR);

        for (Pair<String, ArrayList<String>> sgRNA : sgRNAs) {
            String appendedMotif_ID = sgRNA.getKey();
            ArrayList<String> oligos = sgRNA.getValue();
            System.out.println(appendedMotif_ID + " | " + oligos);
        }
    }

    public static void main(String[] args) throws Exception {
        InvestigatePromoter ip = new InvestigatePromoter();
        ip.initiate();
        ip.run("TERT");
    }
}
