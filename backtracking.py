import argparse, pyBigWig, sys, math, os
from collections import OrderedDict

def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("-r", "--reference", type=str,metavar="reference.bed12")
    base_group.add_argument("-fwd", "--fwd", type=str, nargs='+',metavar="fwd.bw")
    base_group.add_argument("-rev", "--rev", type=str, nargs='+',metavar="rev.bw")
    base_group.add_argument("-o", "--output", type=str)
    base_group.add_argument("-s", "--bgstart", type=int)
    base_group.add_argument("-e", "--bgend", type=int)
    base_group.add_argument("-size","--chromsize",type=str)
    base_group.add_argument("-m", "--multiple", type=int)
    base_group.add_argument("-t", "--time", type=str, nargs='+')
    return parser.parse_args(args)

def strand_division(ref):
    chr_positive = {}
    chr_negative = {}
    chr_positive_sort = {}
    chr_negative_sort = {}
    chrlist = []
    #for i in [x for x in range(1,23)] + ["X","Y"]:
    for i in range(1,6):
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
    for each_key in chr_positive.keys():
        chr_positive_sort[each_key] = OrderedDict(sorted(chr_positive[each_key].items(),key = lambda item:item[1][0]))
    for each_key in chr_negative.keys():
        chr_negative_sort[each_key] = OrderedDict(sorted(chr_negative[each_key].items(),key = lambda item:item[1][0]))
    print chr_positive_sort
    print chr_negative_sort
    return chr_positive_sort,chr_negative_sort

def find_gene_distance(ref):
    """Calculate gene distance"""
    positive,negative = strand_division(ref)
    dFwd = OrderedDict()
    dRev = OrderedDict()
    for each_key in positive.keys():
        lif = OrderedDict(sorted(positive[each_key].items(),key = lambda item:item[1][0]))
        #print(lif)
        length = len(lif.keys())
        #print(length)
        for a in range(length-1):
            distance = int(lif[list(lif.keys())[a+1]][0])-int(lif[list(lif.keys())[a]][0])
            dFwd[list(lif.keys())[a]] = distance
            #print(a)
        dFwd[list(lif.keys())[a+1]] = 100000000
    for each_key in negative.keys():
        lir = OrderedDict(sorted(negative[each_key].items(),key = lambda item:item[1][0]))
        length = len(lir.keys())
        dRev[list(lir.keys())[0]] = 100000000
        for w in range(1,length):
            distance =int(lir[list(lir.keys())[w]][1])-int(lir[list(lir.keys())[w-1]][1])
            dRev[list(lir.keys())[w]] = distance
    return dFwd,dRev

def get_gene_downstream_null(ref,e):
    """Get gene list for no overlapping gene downstream"""
    dFwd,dRev = find_gene_distance(ref)
    lln = []
    for key in dFwd.keys():
        if int(dFwd[key]) > e:
            lln.append(key)
    for key in dRev.keys():
        if int(dRev[key]) > e:
            lln.append(key)
    return lln

def get_open_bw(fw,rw,strand):
    assert strand in "+-"
    sf = []
    sr = []
    if strand == "+":
        for k in range(0,len(fw)):
            hf = 'bw'+str(k)+'forward'
            sf.append(hf)
            sf[k] = pyBigWig.open(fw[k])
    else:
        for k in range(0,len(rw)):
            hr = 'bw'+str(k)+'reverse'
            sr.append(hr)
            sr[k] =pyBigWig.open(rw[k])
    return sf,sr

def fetch_background_signal_downstream_null(fw,rw,s,e,size,chrom,start,end,strand):
    """Fetch GRO-Seq signal in an interval for no overlapping gene downstream"""
    ll_null = []
    chromsizes = {}
    sf,sr = get_open_bw(fw,rw,strand)
    with open(size) as f:
        for line in f:
            l = line.strip("\n").split("\t")
            chromsizes[l[0]] = int(l[1])
    if strand == "+":
        for k in range(0,len(fw)):
            start = start + s
            end = start + e
            if (start+s) > chromsizes[chrom] or (start+e) > chromsizes[chrom]:
                start = chromsizes[chrom] - (e-s)
                end = chromsizes[chrom]
            rd=sf[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            ll_null.append(rd)
    else:
        for k in range(0,len(rw)):
            start = end-e
            end = end-s
            if (end-e)<0 or (end-s)<0:
                start = 0
                end = e-s
            rd = sr[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            ll_null.append(rd)
    return ll_null

def fetch_background_signal_downstream_be(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid):
    """Fetch GRO-Seq signal in an interval for having overlapping gene downstream"""
    ll_be = []
    sf,sr = get_open_bw(fw,rw,strand)
    positive,negative = strand_division(ref)
    if strand == "+":
        for k in range(0,len(fw)):
            idx = list(positive[chrom].keys()).index(Gid)
            start = positive[chrom][list(positive[chrom].keys())[idx+1]][0]-(e-s)
            end = positive[chrom][list(positive[chrom].keys())[idx+1]][0]
            rd = sf[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a] = 0
            ll_be.append(rd)
    else:
        for k in range(0,len(rw)):
            idx = list(negative[chrom].keys()).index(Gid)
            if idx == 0:
                start = 0
                end = e-s
            else:
                start = negative[chrom][list(negative[chrom].keys())[idx-1]][1]
                end = negative[chrom][list(negative[chrom].keys())[idx-1]][1]+(e-s)
                rd = sr[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            ll_be.append(rd)
    return ll_be

def fetch_background_signal(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid):
    """Fetch GRO-Seq signal in an interval"""
    ll = []
    chromsizes = {}
    sf,sr = get_open_bw(fw,rw,strand)
    lln = get_gene_downstream_null(ref,e)
    if Gid in lln:
        ll_null = fetch_background_signal_downstream_null(fw,rw,s,e,size,chrom,start,end,strand)
        ll = ll_null
    else:
        ll_be = fetch_background_signal_downstream_be(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid)
        ll = ll_be
    return ll

def find_mean(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid):
    ll = fetch_background_signal(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid)
    lu = []
    for i in range(len(fw)):
        sum = 0
        for j in range(0,e-s):
            sum = sum + ll[i][j]
        u = sum / (e-s)
        lu.append(u)
    return lu

def find_sd(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid):
    lsd = []
    ll = fetch_background_signal(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid)
    lu = find_mean(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid)
    for i in range(len(fw)):
        sd_element = 0
        for j in range(0,e-s):
            sd_element = sd_element + (ll[i][j]-lu[i])*(ll[i][j]-lu[i])
        sd = math.sqrt(sd_element/(e-s))
        lsd.append(sd)
    return lsd

def fetch_gene_signal(fw,rw,chrom,start,end,strand):
    lsignal = []
    sf,sr = get_open_bw(fw,rw,strand)
    if strand == "+":
        for k in range(0,len(fw)):
            rd=sf[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            lsignal.append(rd)
    else:
        for k in range(0,len(rw)):
            rd = sr[k].values(chrom,start,end)
            for a in range(len(rd)):
                if math.isnan(rd[a]):
                    rd[a]= 0
            rd = rd[::-1]
            lsignal.append(rd)
    return lsignal

def fetch_elongation_signal(fw,rw,chrom,i,window_start,window_end,strand):
    sf,sr = get_open_bw(fw,rw,strand)
    lenlongation = []
    if strand == "+":
        window_start = window_end
        window_end = window_end + 500
        rd=sf[i].values(chrom,window_start,window_end)
        for a in range(len(rd)):
            if math.isnan(rd[a]):
                rd[a]= 0
        for k in range(len(rd)):
            lenlongation.append(rd[k])
    else:
        window_end = window_start
        window_start = window_start - 500
        rd = sr[i].values(chrom,window_start,window_end)
        for a in range(len(rd)):
            if math.isnan(rd[a]):
                rd[a]= 0
        rd = rd[::-1]
        for k in range(len(rd)):
            lenlongation.append(rd[k])
    return lenlongation,window_start,window_end

def find_endpoint(fw,rw,ref,s,e,m,size,chrom,start,end,strand,Gid):
    """Fetch GRO-Seq endpoint relative to TSS"""
    lsignal = fetch_gene_signal(fw,rw,chrom,start,end,strand)
    lu = find_mean(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid)
    lsd = find_sd(fw,rw,ref,s,e,size,chrom,start,end,strand,Gid)
    sf,sr = get_open_bw(fw,rw,strand)
    lendpoint = []
    for i in range(len(fw)):
        window_start = start
        window_end = end
        for j in range(len(lsignal[i])):
            if lsignal[i][j] > (lu[i] + m*lsd[i]):
                break
        for z in range(j,len(lsignal[i])):
            if lsignal[i][z] <= (lu[i] + m*lsd[i]):
                break
        if j == len(lsignal[i])-1:
            endpoint = 0
        elif z == len(lsignal[i])-1:
            a=0
            while True:
                lenlongation,window_start,window_end = fetch_elongation_signal(fw,rw,chrom,i,window_start,window_end,strand)
                for w in range(len(lenlongation)):
                    if lenlongation[w] <= (lu[i] + m*lsd[i]):
                        break
                if w != len(lenlongation)-1:
                    break
                else:
                    a = a+1
                    print(a)
            endpoint = len(lsignal[i])+a*500+w+1
        else:
            endpoint = z
        lendpoint.append(endpoint)
    return lendpoint

def main(args):
    args = parse_args(args)
    ref = args.reference
    fw = args.fwd
    rw = args.rev
    s = args.bgstart
    e = args.bgend
    m = args.multiple
    size = args.chromsize
    o = args.output
    #fl=open('backtrack.txt', 'a')
    fl=open('backtrack2.txt', 'a')
    assert len(fw)==len(rw)
    chrlist = []
    for i in [x for x in range(1,23)] + ["X","Y"]:
        chrlist.append("chr" + str(i))
    with open(ref) as f:
        for line in f:
            l = line.strip("\n").split("\t")
            chrom = l[0]
            if chrom in chrlist:
                start = int(l[1])
                end = int(l[2])
                strand = l[5]
                Gid = l[3]
                lpre=[chrom,l[1],l[2],l[3],l[5]]
                lendpoint = find_endpoint(fw,rw,ref,s,e,m,size,chrom,start,end,strand,Gid)
                for x in lendpoint:
                    lpre.append(x)
                fl.write("\t".join(map(str,lpre))+'\n')
        fl.close()

def run():
    main(sys.argv[1:])

if __name__ == "__main__":
    run()
