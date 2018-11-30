'''
avg_curves.py
@author Laura Colbran, 03/29/2017
some code pulled from EnhancerFinder scripts

given an unlimited set of roc or pr curve files, finds the min, max, and mean
of them, and prints a new curve file.
** assumes all curves have the same set of x-values

USAGE:
python avg_curves.py PATH/TO/CURVES/*curve

'''
import sys

def main():
    files = sys.argv[1:]
    avg = []
    mini = []
    maxi = []
    num = 0
    for f in files:
        curve = open(f, 'r').readlines()
        num += 1
        if len(avg) == 0: #if it's the first file
            avg = curve[2].strip().split(" ")
            for i in range(0,len(avg)):
                avg[i] = float(avg[i])
                mini.append(avg[i])
                maxi.append(avg[i])
        else:
            new_curve = curve[2].strip().split(" ")
            for i in range(0,len(new_curve)):
                y = float(new_curve[i])
                avg[i] += y
                if mini[i] > y: mini[i] = y
                if maxi[i] < y: maxi[i] = y

    # finish average calculation and put into string
    a = ""
    m1 = ""
    m2 = ""
    for i in range(0,len(avg)):
        a += "%.5f " % (avg[i]/num)
        m1 += "%.5f " % mini[i]
        m2 += "%.5f " % maxi[i]

    print "summary_curve\n%s%s\n%s\n%s" % (curve[1],a.strip(),m1.strip(),m2.strip())


if __name__ == "__main__":
  main()
