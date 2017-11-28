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

    public void initiate() throws Exception {
        dgl = new DownstreamGenomicLocus();
        ugl = new UpstreamGenomicLocus();
        tfms = new TFmotifs();
        dgl.initiate();
        ugl.initiate();
        tfms.initiate();
    }

    public void run(String geneName) throws Exception {
        ArrayList<String> upstream = ugl.run(geneName);
        String RightHomologyArm = dgl.run(geneName);
        String LeftHomologyArm = upstream.get(0);
        String PromoterAndFivePrimeUTR = upstream.get(1);
        ArrayList<Pair<String, ArrayList<Integer>>> TranscriptionFactors = tfms.run(PromoterAndFivePrimeUTR);
    }

    public static void main(String[] args) throws Exception {
        InvestigatePromoter ip = new InvestigatePromoter();
        ip.initiate();
        ip.run("TERT");
    }
}
