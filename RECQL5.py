import argparse, pyBigWig, sys, math, os

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-r", "--reference", type=str,metavar="reference.bed12")
    base_group.add_argument("-fwd", "--fwd", type=str, nargs='+',metavar="fwd.bw")
    base_group.add_argument("-rev", "--rev", type=str, nargs='+',metavar="rev.bw")
    base_group.add_argument("-o", "--output", type=str)
    base_group.add_argument("-d", "--depth", type=float)
    base_group.add_argument("-c", "--cons", type=int)
    base_group.add_argument("-t", "--time", type=str, nargs='+')
    base_group.add_argument("-size","--chromsize",type=str)
    return parser.parse_args(args)

def forward_length(fw,d,c,chrom,start,end):
    sf = []
    ll = []
    for k in range(0,len(fw)):
        hf = 'bw'+str(k)+'forward'
        sf.append(hf)
        sf[k] = pyBigWig.open(fw[k]) 
    try:
        for k in range(0,len(fw)):
            rd=sf[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            for i in range(120000-c):
                if max([rd[j] for j in range(i,i+c)]) < d:
                    break
            ll.append(i+start)
        return ll
    except RuntimeError as err:
        pass

def reverse_length(rw,d,c,chrom,start,end):
    sr = []
    ll = []
    for k in range(0,len(rw)):
        hr = 'bw'+str(k)+'reverse'
        sr.append(hr)
        sr[k] = pyBigWig.open(rw[k]) 
    try:
        for k in range(0,len(rw)):
            x=sr[k].values(chrom,end,start)
            rd = x[::-1]
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            for i in range(120000-c):
                if max([rd[j] for j in range(i,i+c)]) < d:
                    break
            ll.append(start-i)
        return ll
    except RuntimeError as err:
        pass

def main(args):
    args = parse_args(args)
    ref = args.reference
    fw = args.fwd
    rw = args.rev
    d = args.depth
    c = args.cons
    t= args.time
    size = args.chromsize
    fl=open('/home/lyy/data1/enlongation/RECQL5/endpoint.txt', 'a')
    chromsizes = {}
    chrlist = []
    for i in [x for x in range(1,23)] + ["X","Y"]:
        chrlist.append("chr" + str(i))
    assert len(fw)==len(rw)==len(t)
    with open(size) as f:
        for line in f:
            l = line.strip("\n").split("\t")
            chromsizes[l[0]] = int(l[1])
    with open(ref) as f:
        for line in f:
            l = line.strip("\n").split("\t")
            if l[5] == "+" :
                chrom = l[0]
                if chrom not in chrlist:
                    pass
                start = int(l[1])
                end = int(l[1])+120000
                if end > chromsizes[chrom]:
                    end = chromsizes[chrom]
                l=[chrom,l[1],l[2],l[3],l[5]]
                for u in forward_length(fw,d,c,chrom,start,end): 
                    l.append(u)
                fl.write("\t".join(map(str,l))+'\n')
            else:
                chrom = l[0]
                if chrom not in chrlist:
                    pass
                start = int(l[2])
                end = int(l[2])-120000
                if end < 0:
                    end = 0
                l=[chrom,l[1],l[2],l[3],l[5]]  
                for u in reverse_length(rw,d,c,chrom,start,end):
                    l.append(u)
                fl.write("\t".join(map(str,l))+'\n')

    fl.close()
    for k in range(0,len(t)):
        sf[k].close()
        sr[k].close()

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()

