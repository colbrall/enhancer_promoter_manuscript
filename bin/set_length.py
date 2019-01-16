''' 
set_length.py
@author Laura Colbran

makes everything the same length by changing the start/end locations
keeps regions centered.

USAGE:
python set_length.py PATH/TO/BED/FILE N

N = desired length (bp)

'''



import sys

def main():
    s = [str.split(line.strip(),'\t') for line in open(sys.argv[1], 'r')]
    target = float(sys.argv[2])/2
    for line in s:
        if line[0].startswith("#"): continue
        cent = float(int(line[2]) - int(line[1]))/2 + int(line[1])
        line[1] = cent - target
        line[2] = cent + target
        print ("%s\t%d\t%d") % (line[0], line[1], line[2])

if __name__== "__main__":
    main()
