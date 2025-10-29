# TOPAS v5.0 cif cleaner
A python script for cleaning unrealistic bond dist. or angle output lines from TOPAS V5. 
Also useful for removing dummy atom sites defined like 'a0' 'a1' 

For bond distance, it simply compares the value calculated by TOPAS and single bond distance calculated based on covalent radii (from wikipedia).
For angle, it removes H centered bond angle values, like O-H-C, as well as any angle value below 90 deg.

NB: alternative workaround is to put 'min_r' and 'max_r' in 'site' declaration line

Python: 3.9+
Dependences: Pandas
