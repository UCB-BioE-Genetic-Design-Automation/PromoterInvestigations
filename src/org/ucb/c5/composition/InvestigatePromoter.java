package org.ucb.c5.composition;

import javafx.util.Pair;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 * @author Manraj Gill
 */
public class InvestigatePromoter {

    private DownstreamGenomicLocus dgl;
    private UpstreamGenomicLocus ugl;
    private TFmotifs tfms;
    private Cas9ConstructOligos c9co;
    private RepairConstructStrategy rcs;
    private AlterationConstructOligos aco;

    public void initiate() throws Exception {
        dgl = new DownstreamGenomicLocus();
        ugl = new UpstreamGenomicLocus();
        tfms = new TFmotifs();
        c9co = new Cas9ConstructOligos();
        rcs = new RepairConstructStrategy();
        aco = new AlterationConstructOligos();
        dgl.initiate();
        ugl.initiate();
        tfms.initiate();
        c9co.initiate();
        rcs.initiate();
        aco.initiate();
    }

    public void run(String geneName, String directory) throws Exception {

        // Run the overarching algorithm
        ArrayList<String> upstream = ugl.run(geneName);
        String RightHomologyArm = dgl.run(geneName);
        String LeftHomologyArm = upstream.get(0);
        String PromoterAndFivePrimeUTR = upstream.get(1);
        ArrayList<ArrayList> TFmotifsInSequence = tfms.run(PromoterAndFivePrimeUTR);
        ArrayList<Pair<String, ArrayList<Integer>>> TranscriptionFactors = TFmotifsInSequence.get(0);
        ArrayList<HashMap> TFinformation = TFmotifsInSequence.get(1);
        ArrayList<Pair<String, ArrayList<String>>> sgRNAs = c9co.run(TranscriptionFactors, PromoterAndFivePrimeUTR);
        HashMap<String, ArrayList<String>> hashMapOfsgRNAs = new HashMap<>();
        for (Pair<String, ArrayList<String>> sgRNA : sgRNAs) {
            String appendedMotif_ID = sgRNA.getKey();
            ArrayList<String> sgRNAoligos = sgRNA.getValue();
            hashMapOfsgRNAs.put(appendedMotif_ID, sgRNAoligos);
        }
        ArrayList<String> repairStrategy = rcs.run(LeftHomologyArm, PromoterAndFivePrimeUTR, RightHomologyArm);
        ArrayList<Pair<String, ArrayList<String>>> alterationOligos = aco.run(TranscriptionFactors, PromoterAndFivePrimeUTR);
        HashMap<String, String> tfNames = TFinformation.get(0);
        HashMap<String, String> tfFamilyNames = TFinformation.get(1);
        HashMap<String, String> consensusSequences = TFinformation.get(2);

        // Parse information from the constructed data structures into a tab separated file
        File file = new File(directory + "/" + geneName);
        FileWriter fw = new FileWriter(file);
        PrintWriter pw = new PrintWriter(fw);
        pw.println("EXPERIMENTAL OVERVIEW:");
        pw.println();
        pw.println("In order to investigate transcription factors in the promoter and 5' UTR of " + geneName + ", use the following two primers to amplify the genomic locus:");
        pw.println();
        pw.println(repairStrategy.get(2) + " and " + repairStrategy.get(3));
        pw.println("The expected size of the amplicon is " + repairStrategy.get(4));
        pw.println("Digest the pCR2.1-TOPO backbone and the amplicon with " + repairStrategy.get(0) + " and " + repairStrategy.get(1) + ", ligate and transform to obtain the repair construct.");
        pw.println();
        pw.println("------------------------------");
        pw.println();
        pw.println("Using this repair construct... For each of the TF motifs identified below, use the Forward_Overlap_Primer and Reverse_Overlap_Primer in an overlab-extension strategy (see outline in README) to devise unique constructs that remove the TF motif.");
        pw.println();
        pw.println("Similarly for each TF motif identified below, use the sgRNA_UP and sgRNA_DOWN in the cloning protocol of px330 constructs.");
        pw.println("The cloning protocol of px330 constructs is adapted from Zhang Lab protocols for CRISPR/Cas9 systems in mammalian cells.");
        pw.println("Detailed protocol available in this package as a PDF title Cas9ConstructProtocol.pdf");
        pw.println();
        pw.println("------------------------------");
        pw.println();
        pw.println("IDENTIFIED TRANSCRIPTION FACTOR BINDING MOTIFS:");
        pw.println();
        pw.println("Count\tMotif_ID\tsgRNA_UP\tsgRNA_DOWN\tForward_Overlap_Primer\tReverse_Overlap_Primer\tLocation\tConsensus_Sequence\tTF_Name\tTF_Family");
        int i = 0;
        for (Pair<String, ArrayList<String>> alteration : alterationOligos) {
            pw.print(i + "\t");

            String appendedMotif_ID = alteration.getKey();
            ArrayList<String> oligos = alteration.getValue();

            // Obtain and print the Motif_ID
            String Motif_ID = appendedMotif_ID.substring(0, 10);
            pw.print(Motif_ID + "\t");

            // Obtain and print the up sgRNA oligos
            ArrayList<String> sgRNAoligos = hashMapOfsgRNAs.get(appendedMotif_ID);
            String sgRNA_UP = sgRNAoligos.get(0);
            pw.print(sgRNA_UP + "\t");
            String sgRNA_DOWN = sgRNAoligos.get(1);
            pw.print(sgRNA_DOWN + "\t");

            // Obtain and print the forward and reverse overlap extension primers
            String Forward_Overlap_Primer = oligos.get(0);
            pw.print(Forward_Overlap_Primer + "\t");
            String Reverse_Overlap_Primer = oligos.get(1);
            pw.print(Reverse_Overlap_Primer + "\t");

            // Obtain and print the Location relative to the Transcriptional Start Site
            String Location = appendedMotif_ID.substring(11);
            int LocationRelativeToTSS = Integer.parseInt(Location) - 1000;
            String LocationToTSS = Integer.toString(LocationRelativeToTSS);
            pw.print(LocationToTSS + "\t");

            // Obtain and print the Consensus Sequence
            String Consensus_Sequence = consensusSequences.get(Motif_ID);
            pw.print(Consensus_Sequence + "\t");

            // Obtain and print the TF Name
            String TF_Name = tfNames.get(Motif_ID);
            pw.print(TF_Name + "\t");

            // Obtain and print the TF Family Name
            String TF_Family = tfFamilyNames.get(Motif_ID);
            pw.print(TF_Family + "\t");

            // Go to the next line for the next motif
            pw.println();
            i++;
        }
        pw.close();
        fw.close();
    }

    public static void main(String[] args) throws Exception {
        Scanner scanner = new Scanner(System.in);
        InvestigatePromoter ip = new InvestigatePromoter();
        ip.initiate();
        for (int i = 0; i < 10; i++) {
            System.out.print("Enter gene name (must have an NCBI gene ID): ");
            String geneName = scanner.next();
            System.out.println("Enter the path of the desired output directory: ");
            String directory = scanner.next();
            ip.run(geneName, directory);
        }
    }
}
