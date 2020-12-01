import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

// example usage 
// java KMers seq/sequences.fa 1 1000
/**
 *
 * @author rgoldst
 */
public class KMers {

    HashMap<String, String> sequenceData = new HashMap<>();
    String refSequence = "";
    int kLength = 15;
    double nSeq = 0.;
    int initialStart = 0;
    int endPos = 1000;
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        KMers kmers = new KMers();
        kmers.readFiles(args);
        kmers.run();
    }
    
    KMers() {
    }

    void readFiles(String[] args){
        if (args.length > 1) {
            initialStart = Integer.parseInt(args[1]);
            endPos = Integer.parseInt(args[2]);
        }
        boolean firstSeq = true;
        try {
            FileReader file = new FileReader(args[0]);
            BufferedReader buff = new BufferedReader(file);
            boolean eof = false;
            while (!eof) {
                String line = buff.readLine();
                if (line == null) {
                    eof = true;
                } else if (firstSeq) {
                    refSequence = buff.readLine();
                    firstSeq = false;
                } else {
                    String name = line.replaceFirst(">", "");
                    String seq = buff.readLine();
                    sequenceData.put(name, seq);
                    System.out.println(name);
                }
            }
        } catch (IOException ioe) {
                System.out.println("Error -- " + ioe.toString());
                System.exit(1);
        }
        nSeq = 1. * sequenceData.size();
    }
    
    void run() {  //the meat
        int nStarts = refSequence.length() - kLength; //length of kmer, 15
        // if this run will encouter the last base of the sequence then
        if(nStarts < endPos) {
            endPos = nStarts;
        }
        for (int iStart = initialStart; iStart < endPos; iStart++) { //take kmer from 1-15, 2-16
            String k = refSequence.substring(iStart, iStart + kLength);
            double nHit = 0.;
            for (String otherSeq : sequenceData.values()) { //look at every other sequence, is that kmer contained in the reference sequence at all?
                if (otherSeq.contains(k)) { // if contains string
                    nHit++;
                } else {
                    boolean intOneMiss = false; // if doesnt contain that string //genius code
                    for (int iSite = 0; iSite < kLength; iSite++) {
                        String oneAway = ".*" + k.substring(0,iSite) + "."; // regex grep type thing
                        if (iSite < kLength-1) {
                            oneAway = oneAway + k.substring(iSite+1, kLength);
                        }
                        oneAway = oneAway + ".*";
                        intOneMiss = intOneMiss || otherSeq.matches(oneAway);
                    } 
                    if (intOneMiss) {
                        nHit++;
                    }
                }
            }
            System.out.println(iStart + "\t" + (nHit/nSeq) + "\t" + k); // what fraction of alt sequences contain that.
        }
    }
}
