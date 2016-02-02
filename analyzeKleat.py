
__version__ = '1.7'

import argparse,os,sys,gzip
from copy import deepcopy
import subprocess
import math

# Additional modules
import pysam

def kleat_int(thing):
    try:
        return int(thing)
    except ValueError:
        return None

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError:
        return '{} already exists. Did not create!'.format(path)

class GTF:
    def __init__(self, seqname=None, source=None, feature=None, start=None, end=None, score=None, strand=None, frame=None, attribute=None):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        try:
            self.start = int(start)
        except (ValueError, TypeError) as e:
            self.start = None
        try:
            self.end = int(end)
        except (ValueError, TypeError) as e:
            self.end = None
        try:
            self.score = int(score)
        except (ValueError, TypeError) as e:
            self.score = None
        self.strand = strand
        self.frame = frame
        self.attribute = attribute

    def __str__(self):
        return ('\t').join([str(x) for x in [self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attribute]])

class BedLine:
    
    def __init__(self, chrom, start, stop, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
        self.chrom = chrom
        self.chromStart = int(start)
        self.chromEnd = int(stop)
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts
        if score:
            self.score = int(score)
        self.strand = strand
        if thickStart:
            self.thickStart = int(thickStart)
        if thickEnd:
            self.thickEnd = int(thickEnd)
        self.itemRgb = itemRgb
        if blockCount:
            self.blockCount = int(blockCount)
        if blockSizes:
            self.blockSizes = [int(x) for x in filter(None, blockSizes.split(','))]
        if blockStarts:
            self.blockStarts = [int(x) for x in filter(None, blockStarts.split(','))]

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand)

class KallistoResult:
    
    def __init__(self, target_id, length, eff_length, est_counts, tpm):
        self.target_id = target_id
        self.length = length
        self.eff_length = eff_length
        self.est_counts = est_counts
        self.tpm = tpm

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}'.format(self.target_id, self.length, self.eff_length, self.est_counts, self.tpm)

class KleatResult:

    def __init__(self, gene, transcript, transcript_strand, coding, contig, chromosome, cleavage_site, within_UTR, distance_from_annotated_site, ESTs, length_of_tail_in_contig, number_of_tail_reads, number_of_bridge_reads, max_bridge_read_tail_length, bridge_read_identities, tail_and_bridge_reads, number_of_link_pairs, max_link_pair_length, link_pair_identities, pas, utr3):
        self.gene = gene
        self.transcript = transcript
        self.transcript_strand = transcript_strand
        self.coding = coding
        self.contig = contig
        self.chromosome = chromosome
        self.cleavage_site = int(cleavage_site)
        if (within_UTR == 'no'):
            self.within_UTR = False
        else:
            self.within_UTR = True
        self.distance_from_annotated_site = kleat_int(distance_from_annotated_site)
        self.ESTs = ESTs
        self.length_of_tail_in_contig = kleat_int(length_of_tail_in_contig)
        self.number_of_tail_reads = kleat_int(number_of_tail_reads)
        self.number_of_bridge_reads = kleat_int(number_of_bridge_reads)
        self.max_bridge_read_tail_length = kleat_int(max_bridge_read_tail_length)
        self.bridge_read_identities = bridge_read_identities
        self.tail_and_bridge_reads = kleat_int(tail_and_bridge_reads)
        self.number_of_link_pairs = kleat_int(number_of_link_pairs)
        self.max_link_pair_length = kleat_int(max_link_pair_length)
        self.link_pair_identities = link_pair_identities
        self.pas = self.utr3 = None
        if (pas != '-'):
            self.pas = [int(x) for x in pas.split(':')]
        if (utr3 != '-'):
            self.utr3 = [int(x) for x in utr3.split('-')]

    def __str__(self):
        atts = [self.gene, self.transcript, self.transcript_strand, self.coding, self.contig, self.chromosome, self.cleavage_site, self.within_UTR, self.distance_from_annotated_site, self.ESTs, self.length_of_tail_in_contig, self.number_of_tail_reads, self.number_of_bridge_reads, self.max_bridge_read_tail_length, self.bridge_read_identities, self.tail_and_bridge_reads, self.number_of_link_pairs, self.max_link_pair_length, self.link_pair_identities, self.pas, self.utr3]
        atts = [str(x) for x in atts]
        return ('\t').join(atts)

def parseConfig(config):
    results = {'kleats': [], 'r1s': [], 'r2s': []}
    header = None
    with open(config, 'r') as f:
        for line in f:
            line = line.strip()
            if line[0] == '[':
                header = line[1:-1]
                if header not in results:
                    sys.exit('Invalid header! Valid headers are: "[kleats]", "[r1s]", and "[r2s]"!')
                continue
            elif not line:
                continue
            results[header].append(line)
    return results

def parse_gene_names(gene_names):
    results = {}
    with open(gene_names,'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            tid, gene = line.strip().split('\t')
            results[tid] = gene
    return results

def parseGTF(gtffile, sources=None, features=None , not_sources=None, not_features=None, gzipped=False, add_chr=True):
    results = []
    if gzipped:
        f = gzip.open(gtffile, 'rb')
    else:
        f = open(gtffile, 'r')
    for line in f:
        if line[0] == '#':
            continue
        gtf = GTF(*line.strip().split('\t'))
        if add_chr:
            gtf.seqname = 'chr' + gtf.seqname
        attributes = {}
        for attr in [x.split() for x in gtf.attribute.split(';')][:-1]:
            attributes[attr[0]] = attr[1][1:-1]
        gtf.attribute = attributes
        if not_sources and gtf.source in not_sources:
            continue
        elif not_features and gtf.feature in not_features:
            continue
        if sources and gtf.source not in sources:
            continue
        elif features and gtf.feature not in features:
            continue
        results.append(gtf)
    f.close()
    return results

def groupGTF(gtf):
    result = {}
    for g in gtf:
        if g.seqname not in result:
            result[g.seqname] = {g.attribute['gene_name']: [g]}
        if g.attribute['gene_name'] not in result[g.seqname]:
            result[g.seqname][g.attribute['gene_name']] = [g]
        else:
            result[g.seqname][g.attribute['gene_name']].append(g)
    for chrom in result:
        for gene in result[chrom]:
            result[chrom][gene].sort(key=lambda x: x.start)
            result[chrom][gene] = mergeGTFList(result[chrom][gene])
    return result

def parseKleat(kleat, min_bridge_read_tail_len=None, min_num_bridge_reads=None, min_tail_len=None, min_num_tail_reads=None, with_pas=False):
    results = []
    with open(kleat, 'r') as f:
        f.readline()
        for line in f:
            result = KleatResult(*line.strip().split('\t'))
            if min_bridge_read_tail_len and (result.max_bridge_read_tail_length < min_bridge_read_tail_len):
                continue
            elif min_num_bridge_reads and (result.number_of_bridge_reads < min_num_bridge_reads):
                continue
            elif with_pas and not result.pas:
                continue
            results.append(result)
    return results

def genTrackLine(name, description=None, _type=None, visibility=2, color=None):
    result = ['name="{}"'.format(name)]
    if description:
        result.append('description="{}"'.format(description))
    if _type:
        result.append('type="{}"'.format(_type))
    if visibility:
        result.append('visibility="{}"'.format(visibility))
    if color:
        result.append('color="{}"'.format(color))
    return 'track ' + (' ').join(result) + '\n'

# Merge gtf line with another gtf line
def mergeTwoGTFs(g1, g2, delim='|'):
    res = GTF()
    res.seqname = g1.seqname
    res.source = g1.source + delim + g2.source
    res.feature = g1.feature + delim + g2.feature
    res.start = min(g1.start, g2.start)
    res.end = max(g1.end, g2.end)
    try:
        res.score = float(g1.score) + g2.score / 2
    except TypeError:
        res.score = None
    res.strand = g1.strand
    res.frame = g1.frame
    res.attribute = g1.attribute
    return res

# list must be sorted
def mergeGTFList(gtf_list):
    res = [gtf_list[0]]
    for i in xrange(1,len(gtf_list)):
        if (gtf_list[i].start <= res[-1].end):
            res[-1] = mergeTwoGTFs(gtf_list[i],res[-1])
        else:
            res.append(gtf_list[i])
    return res

def genFa(clusters, gtf, ref, min_len=20, cap_size=100, delim='|'):
    result = ''
    regions = ''
    intervals = ''
    for chrom in gtf:
        if chrom not in clusters:
            continue
        for gene in gtf[chrom]:
            if gene not in clusters[chrom]:
                continue
            strand = gtf[chrom][gene][0].strand
            cleaved = False
            clip = None
            for i in xrange(len(gtf[chrom][gene])-1):
                d = gtf[chrom][gene][i+1].start - gtf[chrom][gene][i].end
                if d > 3000:
                    clip = i
#                    print gene
#                    print clip
#                    print [[x.start, x.end] for x in gtf[chrom][gene]]
            if clip or clip == 0:
                if strand == '+':
                    gtf[chrom][gene] = gtf[chrom][gene][clip+1:]
                else:
                    gtf[chrom][gene] = gtf[chrom][gene][:clip+1]
                print [[x.start, x.end] for x in gtf[chrom][gene]]
            for region in gtf[chrom][gene]:
                last = region.start
                for kleat in clusters[chrom][gene]:
                    cs = kleat.cleavage_site
                    if cs - last < min_len:
                        continue
                    if last < cs <= region.end:
                        cleaved = True
                        result += (delim).join(['>utr3', gene, strand, chrom, str(last), str(cs)]) + '\n'
                        result += reference.fetch(chrom, last, cs).upper() + '\n'
                        regions += ('\t').join([chrom, str(last), str(cs), gene, '0', strand]) + '\n'
                        intervals += ('\t').join([(delim).join(['utr3', gene, chrom, str(last), str(cs)]), '1', str(cs-last), strand]) + '\n'
                        last = cs+1
                if (region.end - last >= min_len) and cleaved:
                    result += (delim).join(['>utr3', gene, strand, chrom, str(last), str(region.end)]) + '\n'
                    result += reference.fetch(chrom, last, region.end).upper() + '\n'
                    regions += ('\t').join([chrom, str(last), str(region.end), gene, '0', strand]) + '\n'
                    intervals += ('\t').join([(delim).join(['utr3', gene, chrom, str(last), str(region.end)]), '1', str(region.end-last), strand]) + '\n'
            # Capture window after/before end of utr3
            if strand == '+' and cleaved:
                last_exon = gtf[chrom][gene][-1]
                result += (delim).join(['>cap3', gene, strand, last_exon.seqname, str(last_exon.end), str(last_exon.end + cap_size)]) + '\n'
                result += reference.fetch(last_exon.seqname, last_exon.start, last_exon.end + cap_size).upper() + '\n'
                regions += ('\t').join([last_exon.seqname, str(last_exon.end), str(last_exon.end + cap_size), gene, '0', strand]) + '\n'
                intervals += ('\t').join([(delim).join(['cap3', gene, last_exon.seqname, str(last_exon.end), str(last_exon.end + cap_size)]), '1', str(cap_size), strand]) + '\n'
            elif strand == '-' and cleaved:
                last_exon = gtf[chrom][gene][0]
                result += (delim).join(['>cap3', gene, strand, last_exon.seqname, str(last_exon.start - cap_size), str(last_exon.start)]) + '\n'
                result += reference.fetch(last_exon.seqname, last_exon.start - cap_size, last_exon.end).upper() + '\n'
                regions += ('\t').join([last_exon.seqname, str(last_exon.start - cap_size), str(last_exon.start), gene, '0' , strand]) + '\n'
                intervals += ('\t').join([(delim).join(['cap3', gene, last_exon.seqname, str(last_exon.start - cap_size), str(last_exon.start)]), '1', str(cap_size), strand]) + '\n'
    return result.strip(), regions.strip(), intervals.strip()
                
def writeFile(path, name=None, *lines):
    if name:
        result = os.path.join(path,name)
    else:
        result = path
    with open(result, 'w') as f:
        f.write(('\n').join(lines))
    return result

def kallistoIndex(kallisto_path, fasta_path):
    index_path = fasta_path + '.idx'
    index_cmd = [kallisto_path, 'index', '-i', index_path, fasta_path]
    index = subprocess.Popen(index_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return index

def writeKallistoQuant(kallisto_path, index_path, out_path, reads_1, reads_2, bias=False, bootstrap=None, threads=None):
    expression_command = [kallisto_path, 'quant', '-i', index_path, '--pseudobam']
    if bias:
        expression_command.append('--bias')
    if bootstrap:
        expression_command += ['-b', bootstrap]
    if threads:
        expression_command += ['-t', threads]
    expression_command += ['-o', out_path, reads_1, reads_2]
    return (' ').join(expression_command)

def kallistoQuant(kallisto_path, index_path, out_path, reads_1, reads_2, bias=False, bootstrap=None, threads=None, stdout=None):
    expression_command = [kallisto_path, 'quant', '-i', index_path, '--pseudobam']
    if bias:
        expression_command.append('--bias')
    if bootstrap:
        expression_command += ['-b', bootstrap]
    if threads:
        expression_command += ['-t', threads]
    expression_command += ['-o', out_path, reads_1, reads_2]
    if not stdout:
        expression = subprocess.Popen(expression_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        expression = subprocess.Popen(expression_command, stdout=stdout, stderr=subprocess.PIPE)
    return expression

def indexFa(fasta_path):
    try:
        subprocess.call('samtools faidx {}'.format(fasta_path).split())
    except FileNotFoundError:
        return '**WARNING**: Could not locate samtools, unable to index {}'.format(fasta_path)
    return True

def groupKleat(parsed):
    results = {}
    for r in parsed:
        if r.chromosome not in results:
            results[r.chromosome] = {r.gene: [r]}
        if r.gene not in results[r.chromosome]:
            results[r.chromosome][r.gene] = [r]
        else:
            results[r.chromosome][r.gene].append(r)
    return results

def genSponge(gtf, reference, delim='|'):
    result = []
    regions = []
    valid_chroms = ['chr'+str(x) for x in range(1,24)]
    used = set()
    for region in gtf:
        line = ('\t').join([region.seqname, str(region.start), str(region.end), region.attribute['gene_name'], '0', region.strand])
        if (line in used) or (region.seqname not in valid_chroms) or (region.end <= region.start):
            continue
        used.add(line)
        header = (delim).join(['>sponge', region.attribute['gene_name'], region.strand, region.seqname, str(region.start), str(region.end)])
        result.append(header)
        result.append(reference.fetch(region.seqname, region.start, region.end).upper())
        regions.append(line)
    return ('\n').join(result), ('\n').join(regions)

def debug(func):
    def debug_and_call(*args, **kwargs):
        sys.stdout.write()

def cluster(clusters):
    kd = Kdtree(clusters).construct()
    A = kd
    B = min([A.left,A.right], key=lambda x: HAC.centroid(A.loc,x.loc))

def avg(_list):
    _len = len(_list)
    return sum(_list)/_len if _len > 0 else None

def centroid(a,b):
    return abs(avg(a)-avg(b))

def mergeKleatResults(sites):
    d = {'cleavage_site': 0,
         'max_len_br': 0,
         'num_br': 0,
         'len_tail_contig': 0,
         'num_tr': 0}
    count = 0
    for c in sites:
        count += 1
        if c.max_bridge_read_tail_length > d['max_len_br']:
            d['max_len_br'] = c.max_bridge_read_tail_length
        d['num_br'] += c.number_of_bridge_reads
        d['cleavage_site'] += c.cleavage_site
        if c.length_of_tail_in_contig > d['len_tail_contig']:
            d['len_tail_contig'] = c.length_of_tail_in_contig
        d['num_tr'] += c.number_of_tail_reads
    d['cleavage_site'] /= count
    res = deepcopy(sites[0])
    res.cleavage_site = d['cleavage_site']
    res.length_of_tail_in_contig = d['len_tail_contig']
    res.number_of_bridge_reads = d['num_br']
    res.number_of_tail_reads = d['num_tr']
    res.max_bridge_read_tail_length = d['max_len_br']
    return res

def linkage(sites, window=20):
    length = len(sites)
    if length > 1:
        _min = float('inf')
        r = s = None
        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                dist = abs(sites[i] - sites[j])
                if dist < _min:
                    r = i
                    s = j
                    _min = dist
        if _min > window:
            return sites
        sites[r] = (sites[r] + sites[s])/2
        del(sites[s])
        linkage(sites, _min)
    return sites

def kleatLinkage(sites, window=20):
    length = len(sites)
    if length > 1:
        _min = float('inf')
        r = s = None
        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                dist = abs(sites[i].cleavage_site - sites[j].cleavage_site)
                if dist < _min:
                    r = i
                    s = j
                    _min = dist
        #sites[r] = _min
        if _min <= window:
            sites[r] = mergeKleatResults([sites[r],sites[s]])
            del(sites[s])
            kleatLinkage(sites)
    return sites

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Take results from KLEAT and generate fasta sequences for 3\'UTR isoforms')
    
    parser.add_argument('config', help='A configuration file containing paths to kleat outputs, and reads for each dataset')
    parser.add_argument('-features', default='/projects/dmacmillanprj2/polya/encode/analysis/fasta_generation/analyzeKleat/Homo_sapiens.GRCh37.75.gtf', help='A file containing genomic features in GTF format (Can be downloaded from Ensembl)')
    parser.add_argument('-gene_names', default='/home/dmacmillan/annotations/ensembl/ensemblToGeneName.original', help='A two column TSV file with annotation labels on the left and corresponding gene names on the right. If header is present, it must start with a \# sign')
    parser.add_argument('-reference', default='/home/dmacmillan/references/hg19/hg19.fa', help='Reference genome in fasta format (index must be in same folder). Default is \'/home/dmacmillan/references/hg19/hg19.fa\'')
    parser.add_argument('-kallisto', default='/gsc/btl/linuxbrew/bin/kallisto', help='Path to Kallisto executable. Default is \'/gsc/btl/linuxbrew/bin/kallisto\'')
    parser.add_argument('-de', '--delim', default='|', help='Set the delimiter for headers in fasta files. This should be a character that is not contained within any gene names, transcript id\'s, or other names. Default is \'|\' neglecting surrounding quotations')
    parser.add_argument('-cw', '--cluster_window', type=int, default=20, help='Set the window size for clustering KLEAT cleavage sites. Default = 20')
    parser.add_argument('-d', '--debug', action='store_true', help='Print detailed debugging information')
    parser.add_argument('-m', '--min_seq_len', type=int, default=20, help='The minimum length of sequence to be output in the fasta file. Default is 20')
    #parser.add_argument('-n', '--name', default=None, help='Name for dataset. Default will be the basename of all -ks arguments delimited by the --delim parameter.')
    parser.add_argument('-bi', '--bias', action='store_false', help='Bias parameter for Kallisto, enabled by default. Use this flag to disable')
    parser.add_argument('-bo', '--bootstrap', default='100', help='Bootstrap parameter for Kallisto. Default = 100')
    parser.add_argument('-t', '--threads', default='8', help='Number of threads to use for Kallisto. Default = 8')
    parser.add_argument('-c', '--cap_size', type=int, default=100, help='Length of region 3\'-adjacent to 3\'UTR to capture. Default = 100')
    parser.add_argument('-mbrtl', '--min_bridge_read_tail_len', type=int, default=2, help='Disregard KLEAT calls which have a maximum bridge read tail length less than this number. Default is 2')
    parser.add_argument('-mnbr', '--min_num_bridge_reads', type=int, default=2, help='Disregard KLEAT calls which have a number of bridge reads less than this number. Default is 2')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output to. Default is current directory')
    parser.add_argument('-nq', '--no_quant', action='store_true', help='Set this flag to skip Kallisto quant phase')

    args = parser.parse_args()

    config = parseConfig(args.config)

    mass = []

    for dataset in config['kleats']:
        if args.debug:
            sys.stdout.write('Parsing kleat {}...'.format(dataset))
        parsed = parseKleat(dataset, min_num_bridge_reads=args.min_num_bridge_reads, min_bridge_read_tail_len=args.min_bridge_read_tail_len, with_pas=True)
        mass += parsed
        if args.debug:
            print 'DONE'

    # Sort
    if args.debug:
        sys.stdout.write('sorting kleat results...')
    mass = sorted(mass, key=lambda x: x.cleavage_site)
    if args.debug:
        print 'DONE'

    # Group
    if args.debug:
        sys.stdout.write('grouping kleat results...')
    grouped = groupKleat(mass)
    if args.debug:
        print 'DONE'
        
    for chrom in grouped:
        for gene in grouped[chrom]:
            sites = grouped[chrom][gene]
            if args.debug:
                print 'Clustering {} {} {}'.format(chrom, gene, [x.cleavage_site for x in sites])
            sites = kleatLinkage(sites, args.cluster_window)
            if args.debug:
                print [x.cleavage_site for x in sites]

    # Global stuff
    # Parse gene names
    if args.debug:
        sys.stdout.write('Parsing gene names...')
    gene_names = parse_gene_names(args.gene_names)
    if args.debug:
        print 'DONE'
    if args.debug:
        sys.stdout.write('Parsing 3\'UTR sequences...')
    gtf = parseGTF(args.features, sources=['protein_coding'], features='UTR')
    if args.debug:
        print 'DONE'
    if args.debug:
        sys.stdout.write('Parsing sponge sequences...')
    #sponges = parseGTF(args.features, sources=['protein_coding'], features='CDS')
    sponges = parseGTF(args.features, not_sources=['protein_coding'], not_features='UTR')
    if args.debug:
        print 'DONE'
    if args.debug:
        sys.stdout.write('Grouping 3\'UTRs...')
    ggtf = groupGTF(gtf)
    if args.debug:
        print 'DONE'
    if args.debug:
        sys.stdout.write('Loading reference...')
    reference = pysam.FastaFile(args.reference)
    if args.debug:
        print 'DONE'
    if args.debug:
        sys.stdout.write('Generating 3\'UTR regions...')
    utr3_fasta, utr3_regions, utr3_intervals = genFa(grouped, ggtf, reference, cap_size=args.cap_size, delim=args.delim)
    utr3_regions = genTrackLine('3\'UTRs', description='3\'UTR sequences', color='50,50,255') + utr3_regions
    if args.debug:
        print 'DONE'
    if args.debug:
        sys.stdout.write('Generating sponge sequences...')
    sponge_fasta, sponge_regions = genSponge(sponges, reference)
    sponge_regions = genTrackLine('Sponge', description='non 3\'utr sequences', color='0,100,0') + sponge_regions
    if args.debug:
        print 'DONE'
    
    utr3_intervals_path = os.path.join(args.outdir, 'utr3_intervals')
    sponge_regions_path = os.path.join(args.outdir, 'sponge_regions.bed')
    utr3_fasta_path = os.path.join(args.outdir, 'utr3_regions.fa')
    utr3_regions_path = os.path.join(args.outdir, 'utr3_regions.bed')
    writeFile(utr3_fasta_path, None, utr3_fasta, sponge_fasta)
    writeFile(utr3_regions_path, None, utr3_regions)
    writeFile(sponge_regions_path, None, sponge_regions)
    writeFile(utr3_intervals_path, None, utr3_intervals)
    if args.debug:
        print '3\'UTR fasta file written to: {}'.format(utr3_fasta_path)
        print '3\'UTR regions file written to: {}'.format(utr3_regions_path)
        print 'Sponge regions file written to: {}'.format(sponge_regions_path)

    if args.debug:
        print 'samtools faidx {}...'.format(utr3_fasta_path)
    res = indexFa(utr3_fasta_path)
    print res

    index_path = utr3_fasta_path + '.idx'

    if not os.path.isfile(index_path):
        index = kallistoIndex(args.kallisto, utr3_fasta_path)
        if args.debug:
            sys.stdout.write('Indexing {}...'.format(utr3_fasta_path))
        o,e = index.communicate()
        if args.debug:
            print 'DONE'
            print 'Index written to: {}'.format(index_path)
            writeFile(args.outdir, 'kallisto.index.o', o)
            writeFile(args.outdir, 'kallisto.index.e', e)

    if args.no_quant:
        sys.exit()

    for i in xrange(len(config['r1s'])):
        sample = config['kleats'][i].split('.')[0]
        sam_path = os.path.join(args.outdir, 'alignment.sam')
        quant = kallistoQuant(args.kallisto, index_path, args.outdir, config['r1s'][i], config['r2s'][i], bias=args.bias, bootstrap=args.bootstrap, threads=args.threads)
        if args.debug:
            print 'Quantifying...'
        quant = quant.communicate()
        if args.debug:
            print 'DONE'
        tsv_path = os.path.join(args.outdir, 'abundance.tsv')
        rename = os.path.join(args.outdir, os.path.basename(config['kleats'][i]).split('.')[0] + '.tsv')
        os.rename(tsv_path, rename)
