import simplejson

file = open('../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf')
output = open('CHBJPTlowcovhash','w')
hash = {}
lines = file.readlines(1000000)
while(lines!=[]):
    for l in lines:
        if not l.startswith('#'):
            tokens = l.strip('\n').split('\t')
            key = 'chr'+tokens[0]+'pos'+tokens[1]
            hash[key] = []
            for t in tokens[9:]:
                if t!= ".":
                    matrixentry = int(t[0]) + int(t[2])
                    hash[key].append(matrixentry)
    lines = file.readlines(1000000)

simplejson.dump(hash, output)
