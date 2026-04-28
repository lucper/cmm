# How to run SLIDER

Unzip the `SliderLight.zip` file and `cd` to the `SliderLight` directory. Then, run:
```bash
java -cp dist/SliderLight.jar Framework.Framework
```
This command should yield:
```text
Arguments l, d, seq, int, o and m are required
Example usage: java Framework -l 8 -d 3 -seq seq.fasta -int prot.int -o output -m HCThreaded

Required options:
-l      The length of the desired motifs
-d      The amount of wildcards in the desired motifs
-seq    The fasta file containing the sequence information
-int    The file containing the interaction information
-o      The output file
-m      The desired method (HCTimed, BruteForceThreaded,slider++)

Other options:
-t      The amount of threads to use (standard uses as many as there are processors available
-a      The amount of best scores to keep (standard 1,000)
-f      The neighbor function to use for HCTimed
-st     allows for different scoretypes (binomial, cover, difference, p, v, weightedv, x2)
-min    the amount of minutes to run (used only for timed methods)
```
