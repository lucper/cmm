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

# How to run D-STAR

Unzip the `DSTAR.zip` file and `cd` to the `Program` directory.
The executables therein are 32-bit x86 binaries, so we need to install the following dependency.
For Ubuntu/Debian systems, run:
```bash
sudo dpkg --add-architecture i386
sudo apt-get update
sudo apt-get install libc6-i386
```
Then test it:
```bash
chmod +x DSTAR.LinuxRedHat
./DSTAR.LinuxRedHat
```
The secon command should yield:
```text
File not found - (null)
File not found - SHELL=/bin/bash
Current Usage = 16040 bytes, Max = 36056
Total Time 0.000174 s
```
(The numbers may be different in your machine.)
