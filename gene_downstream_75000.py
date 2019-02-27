import argparse, pyBigWig, sys, math, os
from collections import OrderedDict

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-r", "--reference", type=str,metavar="reference.bed12")
    return parser.parse_args(args)

def strand_division(ref):
    chr_positive = {}
    chr_negative = {}
    chrlist = []
    for i in [x for x in range(1,23)] + ["X","Y"]:
        chrlist.append("chr" + str(i))
        chr_positive["chr" + str(i)] = {}
        chr_negative["chr" + str(i)] = {}
    with open(ref,"r ") as file:
        for each_line in file:
            line = each_line.rstrip().split("\t")
            if line[0] not in chrlist:
                pass
            elif line[5] == "+":
                chr_positive[line[0]][line[3]] = [int(line[1]),int(line[2])]
            elif line[5] == "-":
                chr_negative[line[0]][line[3]] = [int(line[1]),int(line[2])]
    return chr_positive,chr_negative

def find_gene_distance(ref):
    positive,negative = strand_division(ref)
    dFwd = OrderedDict()
    dRev = OrderedDict()
    for each_key in positive.keys():
        lif = OrderedDict(sorted(positive[each_key].items(),key = lambda item:item[1][0]))
        length = len(lif.keys())
        for i in range(length-1):
            distance = int(lif[list(lif.keys())[i+1]][0])-int(lif[list(lif.keys())[i]][0])
            dFwd[list(lif.keys())[i]] = distance
        dFwd[list(lif.keys())[i+1]] = 100000000
    for each_key in negative.keys():
        lir = OrderedDict(sorted(negative[each_key].items(),key = lambda item:item[1][0]))
        length = len(lir.keys())
        dRev[list(lir.keys())[0]] = 100000000
        for i in range(1,length):
            distance =int(lir[list(lir.keys())[i]][1])-int(lir[list(lir.keys())[i-1]][1])
            dRev[list(lir.keys())[i]] = distance
    return dFwd,dRev

def main(args):
    args = parse_args(args)
    ref = args.reference
    fl=open('gene_75000.txt', 'a')
    dFwd,dRev = find_gene_distance(ref)
    ll = []
    for key in dFwd.keys():
        if int(dFwd[key]) > 75000:
            ll.append(key)
    for key in dRev.keys():
        if int(dRev[key]) > 75000:
            ll.append(key)
    print len(ll)
    with open(ref,"r ") as file:
        for each_line in file:
            line = each_line.rstrip().split("\t")
            if line[3] in ll:
                fl.write(each_line)
    fl.close()

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
