package org.ucb.c5.composition;

import javafx.util.Pair;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.io.File;
import org.ucb.c5.utils.FileUtils;
import java.util.regex.*;


/**
 * Data of Homo Sapiens Transcription Factor motifs from the CIS-BP database:
 *
 * "Determination and inference of eukaryotic transcription factor sequence specificity."
 * Weirauch MT, Yang A, Albu M, Cote AG, Montenegro-Montero A, Drewe P, Najafabadi HS,
 * Lambert SA, Mann I, Cook K, Zheng H, Goity A, van Bakel H, Lozano JC, Galli M, Lewsey
 * MG, Huang E, Mukherjee T, Chen X, Reece-Hoyes JS, Govindarajan S, Shaulsky G, Walhout
 * AJ, Bouget FY, Ratsch G, Larrondo LF, Ecker JR, Hughes TR.
 * Cell. 2014 Sep 11;158(6):1431-43. doi: 10.1016/j.cell.2014.08.009.
 * PMID: 25215497
 * URL: http://cisbp.ccbr.utoronto.ca/index.php
 *
 * Initiation of org.ucb.c5.composition.TFmotifs parses the data from CIS-BP database
 * Running of org.ucb.c5.composition.TFmotifs identifies TF binding sites in the input
 * promoter and 5' UTR sequence
 *
 * @author Manraj Gill
 */

public class TFmotifs {

    private HashMap<String, String> tfIDs;
    private HashMap<String, String> tfNames;
    private HashMap<String, String> dbIDs;
    private HashMap<String, String> consensusSequences;

    public void initiate() throws Exception {
        // Read the TF_Information_all_motifs.txt file
        //    [0] TF_ID: internal CIS-BP ID for the TF
        //    [3] Motif_ID: internal CIS-BP ID for associated motif
        //    [9] TF_Name: name of TF family
        //    [13] DBID: motif ID from associated study
        // Store these options as a key-value (K,V) paired HashMaps named tfIDs, tfNames and dbIDs
        //    K: Unique Motif_ID
        //    V: Corresponding TF_ID, TF_Name and DBID
        String TF_information = FileUtils.readResourceFile("Homo Sapiens TF Motifs from CIS-BP Database/TF_Information_all_motifs_plus.txt");
        String[] rows = TF_information.split("\\r|\\r?\\n");
        tfIDs = new HashMap<>();
        tfNames = new HashMap<>();
        dbIDs = new HashMap<>();
        for (int i = 0; i < rows.length; i++) {
            String row = rows[i];
            String[] columns = row.split("\t");
            tfIDs.put(columns[3], columns[0]);
            tfNames.put(columns[3], columns[9]);
            dbIDs.put(columns[3], columns[13]);
        }

        // Read the PositionWeightMatrices (PWMs) files, each file is named by its Motif_ID.
        // For each file, create a HashMap with the name of the file (the Motif_ID) as the key
        // and the parsed position weight matrix as a string. The parsing of the position weight
        // matrices is based on the frequencies of the nucleotide bases found at each position.
        consensusSequences = new HashMap<>();
        String pwmsPath = "/Users/Manraj/Documents/GitHub/PromoterInvestigations/src/org/ucb/c5/Homo Sapiens TF Motifs from CIS-BP Database/PositionWeightMatrices";
        File pwmsDirectory = new File(pwmsPath);
        File[] pwms = pwmsDirectory.listFiles();
        for (File pwm : pwms) {
            if (pwm.isFile()) {
                // Obtain the path of the file with the Position Weight Matrix
                // and save it locally along with its Motif_ID (w/o .txt)
                String path = "Homo Sapiens TF Motifs from CIS-BP Database/PositionWeightMatrices/" + pwm.getName();
                String Motif_ID = pwm.getName().substring(0, 10);
                // Create a local string called consensusSequence that corresponds
                // to the consensus sequence for the associated Motif_ID
                String consensusSequence = "";
                // Read the file's contents and create the consensusSequence string
                // based on the frequencies of the nucleotides. If no base is found
                // at a frequency greater than 0.5, add 'n' to the consensus sequence.
                // If the position weight matrix is not populated for this file, then
                // ignore and move on without adding it to the HashMap.
                String motif = FileUtils.readResourceFile(path);
                String[] positions = motif.split("\\r|\\r?\\n");
                if (positions.length == 1) {
                    continue;
                }
                int countSpecific = 0;
                for (int i = 1; i < positions.length; i++) {
                    String position = positions[i];
                    String[] bases = position.split("\t");
                    String baseAtThisPosition = null;
                    double Afrequency = Double.parseDouble(bases[1]);
                    if (Afrequency > 0.5) {
                        baseAtThisPosition = "A";
                        countSpecific += 1;
                    }
                    double Cfrequency = Double.parseDouble(bases[2]);
                    if (Cfrequency > 0.5) {
                        baseAtThisPosition = "C";
                        countSpecific += 1;
                    }
                    double Gfrequency = Double.parseDouble(bases[3]);
                    if (Gfrequency > 0.5) {
                        baseAtThisPosition = "G";
                        countSpecific += 1;
                    }
                    double Tfrequency = Double.parseDouble(bases[4]);
                    if (Tfrequency > 0.5) {
                        baseAtThisPosition = "T";
                        countSpecific += 1;
                    }
                    if (baseAtThisPosition == null) {
                        baseAtThisPosition = ".";
                    }
                    consensusSequence += baseAtThisPosition;
                }
                if ((countSpecific * 100 / positions.length) < 50) {
                    continue;
                }
                consensusSequences.put(Motif_ID, consensusSequence);
            }
        }
    }

    public ArrayList<Pair<String, ArrayList<Integer>>> run(String promoterAndFivePrimeUTR) throws Exception {
        ArrayList<Pair<String, ArrayList<Integer>>> TFmotifsInSequence = new ArrayList<>();
        for (HashMap.Entry<String, String> entry : consensusSequences.entrySet()) {
            String Motif_ID = entry.getKey();
            String consensusSequence = entry.getValue();
            Pattern pattern = Pattern.compile(consensusSequence);
            Matcher matcher = pattern.matcher(promoterAndFivePrimeUTR);
            ArrayList<Integer> matchStartIndices = new ArrayList<>();
            while (matcher.find()) {
                int indexOfMatch = matcher.start();
                matchStartIndices.add(indexOfMatch);
                System.out.println(tfNames.get(Motif_ID) + " | " + indexOfMatch);
            }
            Pair<String, ArrayList<Integer>> matches = new Pair<>(Motif_ID, matchStartIndices);
            TFmotifsInSequence.add(matches);
        }
        return TFmotifsInSequence;
    }

    public static void main(String[] args) throws Exception {
        TFmotifs tfm = new TFmotifs();
        tfm.initiate();
        // TERT's upstream genomic locus
        ArrayList<Pair<String, ArrayList<Integer>>> TFmotifsInTERT = tfm.run("AGTGGCCGTGTGGCTTCTACTGCTGGGCTGGAAGTCGGGCCTCCTAGCTCTGCAGTCCGAGGCTTGGAGCCAGGTGCCTGGACCCCGAGGTTGCCCTCCACCCTGTGCGGGCGGGATGTGACCAGATGTTGGCCTCATCTGCCAGACAGAGTGCCGGGGCCCAGGGTCAAGGCCGTTGTGGCTGGTGTGAGGCGCCCGGTGCGCGGCCAGCAGGAGCGCCTGGCTCCATTTCCCACCCTTTCTCGACGGGACCGCCCCGGTGGGTGATTAACAGATTTGGGGTGGTTTGCTCATGGTGGGGACCCCTCGCCGCCTGAGAACCTGCAAAGAGAAATGACGGGCCTGTGTCAAGGAGCCCAAGTCGCGGGGAAGTGTTGCAGGGAGGCACTCCGGGAGGTCCCGCGTGCCCGTCCAGGGAGCAATGCGTCCTCGGGTTCGTCCCCAGCCGCGTCTACGCGCCTCCGTCCTCCCCTTCACGTCCGGCATTCGTGGTGCCCGGAGCCCGACGCCCCGCGTCCGGACCTGGAGGCAGCCCTGGGTCTCCGGATCAGGCCAGCGGCCAAAGGGTCGCCGCACGCACCTGTTCCCAGGGCCTCCACATCATGGCCCCTCCCTCGGGTTACCCCACAGCCTAGGCCGATTCGACCTCTCTCCGCTGGGGCCCTCGCTGGCGTCCCTGCACCCTGGGAGCGCGAGCGGCGCGCGGGCGGGGAAGCGCGGCCCAGACCCCCGGGTCCGCCCGGAGCAGCTGCGCTGTCGGGGCCAGGCCGGGCTCCCAGTGGATTCGCGGGCACAGACGCCCAGGACCGCGCTTCCCACGTGGCGGAGGGACTGGGGACCCGGGCACCCGTCCTGCCCCTTCACCTTCCAGCTCCGCCTCCTCCGCGCGGACCCCGCCCCGTCCCGACCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCCCTTCCTTTCCGCGGCCCCGCCCTCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGCTGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCCCGCG");
    }
}
