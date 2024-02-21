import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class seqkitScript {
    public static void main(String[] args) {
        // bedtools method
        List<String> list = new ArrayList<String>();
        list.add("awk -F'\\\\t' '$3 == \\\"repeat_region\\\" {split($9, a, \\\";\\\"); for (i in a) {split(a[i], b, \\\"=\\\"); if (b[1] == \\\"ID\\\") TE_id = b[2];}  print $1, $4-1, $5, TE_id, \\\".\\\", $7}' chr1A.gff3 > annotations.bed");
        list.add("gff2bed <JM44.repeat.masked.gff > annotation.bed");
        list.add("bedtools getfasta -s -fi JM44.repeat.masked.fasta -bed annotation.bed -fo teseq.fasta");
        list.add("seqkit stats teseq.fasta");
        list.add("python DataProcessing.py > call.txt");
        list.add("bedtools merge -i annotation.bed > merge.bed");
        String[] seqProcessing = list.toArray(new String[list.size()]);
        runCommands(seqProcessing);
        // seqkit to parent
        List<String> list1 = new ArrayList<>();
        list1.add("cat annotation.bed | grep -w 'Parent' > parent.bed");
        list1.add("bedtools getfasta -s -fi 01data/JM44.repeat.masked.fasta -bed new.bed -fo parent.fasta");
        list1.add("cat annotation.bed | grep -v 'Parent' > vParent.bed");
        list1.add("bedtools getfasta -s -fi 01data/JM44.repeat.masked.fasta -bed new.bed -fo vParent.fasta");
        list1.add("seqkit stats 02result/teSeq.fasta parent.fasta vParent.fasta");
        list1.add("cat 01data/");
        list1.add("seqkit fx2tab -nliH parent.fasta -o new.txt");
        list1.add("python DataProcessing.py > table.txt");
        String[] parent = list1.toArray(new String[list1.size()]);
        runCommands(seqProcessing);
        // seq analysis of structural
        List<String> list2 = new ArrayList<>();
        list2.add("cat clrh.bed | grep -w \"Copia_LTR_retrotransposon\" | wc -l");
        String[] structural = list2.toArray(new String[list2.size()]);
    }

    public static void runCommands(String[] commands) {
        try {
            // build ProcessBuilder
            ProcessBuilder processBuilder = new ProcessBuilder();
            for (String command : commands) {
                processBuilder.command(command.split(" "));
            }
            // boot process
            Process process = processBuilder.start();
            // get inputStream of process
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
            // read the output of the command
            String line;
            while ((line = reader.readLine()) != null) {
                System.out.println(line);
            }
            // wait for the finish of process
            int exitCode = process.waitFor();
            System.out.println("Exit Code: " + exitCode);
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }
}