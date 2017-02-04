import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
import pandas as pd
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from StringIO import StringIO
import time

rase_root_dir = '/home/miladim/repositories/RaSE/'
rase_src_dir = os.path.join(rase_root_dir, 'code')
sys.path = [rase_src_dir] + sys.path

eden_root_dir = '/home/miladim/repositories/EDeN/'
eden_src_dir = os.path.join(eden_root_dir)
sys.path = [eden_src_dir] + sys.path

from RaSE import make_fold, make_fold_vectorize

def main(argv):

    input_file = sys.argv[1]
    output_file_prefix = sys.argv[2]
    
    start_time = time.time()

    #initiate empty dataframes for RNAsnp
    df_rnasnp= pd.DataFrame(columns=['SNP','w','Slen','GC','interval','d_max','p-value','interval.1','r_min','p-value.1','ID'])

    #initiate empty dataframes for remuRNA
    df_remurna= pd.DataFrame(columns=['SNP','MFE(wt)','MFE(mu)','dMFE','H(wt||mu)','GCratio','ID'])
    
    rase_scores =[]
    lcount = 0
    fasta_sequences = SeqIO.parse(open(input_file),'fasta')
    for fasta in fasta_sequences:
        lcount += 1
        print '\r{}..' .format(lcount), 
        id, desc, sequence = fasta.id,fasta.description, str(fasta.seq)
        #extract snp from description or id
        snp = [desc[len(desc)-5] + "201"+ desc[len(desc)-1]]
        print snp,
        tmp_seq_fa = NamedTemporaryFile(suffix='.fa', delete=False)
        tmp_seq_fa.write(">" +desc + "\n")
        tmp_seq_fa.write(sequence)
        tmp_seq_fa.close()

        #run RNAsnp
        res_rnasnp=run_RNAsnp(tmp_seq_fa.name,snp,200)
        res_rnasnp=res_rnasnp.assign(ID=id)
        df_rnasnp = df_rnasnp.append(res_rnasnp, ignore_index=True)
        
        #run remuRNA
        res_remurna=run_remuRNA(tmp_seq_fa.name,snp)
        res_remurna=res_remurna.assign(ID=id)
        df_remurna = df_remurna.append(res_remurna, ignore_index=True)

        ##run RaSE
        rase_scores=rase_scores+ [[id,snp[0],run_RaSE(sequence,snp, window=200)]]
                    
        # remove temp file
        os.remove(tmp_seq_fa.name)
        
    df1 = df_rnasnp.set_index('ID')
    df1['tool-parameters:windows|size=200']=''   
    df1.to_csv(path_or_buf=output_file_prefix+"_rnasnp.csv",sep="\t")
    
    df2 = df_remurna.set_index('ID')
    df2['tool-parameters:']=''
    df2.to_csv(path_or_buf=output_file_prefix+"_remurna.csv",sep="\t")
    
    df_rase= pd.DataFrame(rase_scores,columns=['ID','SNP','Score'])
    df_rase =df_rase.set_index('ID')
    df_rase['tool-parameters:window=200|avg_bp_prob_cutoff=0.01|hard_threshold=0.5|max_num_edges=3']='' 
    df_rase.to_csv(path_or_buf=output_file_prefix+"_rase.csv",sep="\t")
    

    print("--- %s seconds ---" % (time.time() - start_time))

def run_remuRNA(wild_fa, snp_tags, window=None):
    """
    A python wrapper invoking remuRNA tool.
    Please check remuRNA documentation for details.
    Call example: run_remuRNA('./wild.fa', ['G20C'])
    Parameters
    ----------
    wild_fa : str
        Fasta file name containing one RNA sequence
    snp_tags : list
        Set of SNP tags required to be evaluatued on the input sequence.
        Warning: remuRNA accepts only a single tag in each call.

    Returns
    -------
    dataframe
        Pandas table of standard output

    """
    assert(len(snp_tags)==1)
    if not os.path.isfile(wild_fa):
        raise RuntimeError ("Input fasta %s does not exist" % in_fasta_file)

    # Write RNA sequence and SNP tags to a temporary file, TODO: Remove the temporary file?
    tmp_fa = NamedTemporaryFile(suffix='.fa', delete=False)
    with open (wild_fa) as in_fa_handle:
        for line in in_fa_handle:
            tmp_fa.write(line)
    tmp_fa.write('\n'.join(['*'+tag for tag in snp_tags]))

    tmp_fa.close()

    # Make a shell command line
    cmd = '$(which remuRNA) {} -p=4 '.format(tmp_fa.name)
    if window is not None:
        cmd += '-w={}'.format(int(window))
    # print cmd
    p = Popen( cmd , stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if err:
        raise RuntimeError("Error in calling remuRNA\n{}\n{}\n".format(out, err))


    os.remove(tmp_fa.name)

    # print out
    df = pd.read_table(StringIO(out))
    return  df
    #return out

def run_RNAsnp(wild_fa, snp_tags, window=None):
    """
    A python wrapper invoking RNAsnp tool.
    Please check RNAsnp documentation for details.

    Call example: run_RNAsnp('./wild.fa', ['G20C'])
    Parameters
    ----------
    wild_fa : str
        Fasta file name containing one RNA sequence
    snp_tags : list
        Set of SNP tags required to be evaluatued on the input sequence
    window : int
        If None, the RNAsnp window (-w) size. Windows larger than 800 will be passed as 800.

    Returns
    -------
    dataframe
        Pandas table of standard output

    """

    # Write SNP tags to a temporary file, TODO: Remove the temporary file?
    snp_file = NamedTemporaryFile(delete=False)
    snp_file.write('\n'.join(snp_tags))
    snp_file.close()

    if not os.path.isfile(wild_fa):
        raise RuntimeError ("Input fasta %s does not exist" % in_fasta_file)

    # Make a shell command line
    cmd = 'RNAsnp -f {} -s {} '.format(wild_fa, snp_file.name)
    if window is not None:
        if window > 800:
            print "WARNING RNAsnp window reduced to max possible: 800"
            window = 800
        cmd += '-w {}'.format(int(window))
    # print cmd
    p = Popen( cmd , stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if err:
        #raise RuntimeError("Error in calling RNAsnp\n{}\n{}\n".format(out, err))
        print "Error in calling RNAsnp\n{}\n{}\n".format(out, err)
    # print out
    
    os.remove(snp_file.name)
        
    out_cleaned = ""
    for line in out.split('\n'):
        if 'error' in line.lower():
            raise RuntimeError("RNASNP returned error for: {} message is:{}".format(wild_fa, line))
        elif 'warning' in line.lower():
            print ("ERROR: RNASNP returned warning for: {} message is:{}".format(wild_fa, line))
        else:
            out_cleaned += line+"\n"
    
    df_RNAsnp = pd.read_table(StringIO(out_cleaned))
    return  df_RNAsnp #.add_suffix('RNAsnp:')
    #return out_cleaned


def run_RaSE(wild_seq, snp_tags, window=150, avg_bp_prob_cutoff=0.01,
                          hard_threshold=0.5,
                          max_num_edges=3,):
    assert(len(snp_tags)==1)
    import re
    
    matches = re.match(r'(\D)(\d+)(\D)', snp_tags[0])
    tag_tup = matches.groups()
    tag_tup = (tag_tup[0], int(tag_tup[1])-1, tag_tup[2]) 
    
    fold = make_fold(window_size=window,
                          max_bp_span=int(window*(0.80)),
                          avg_bp_prob_cutoff=avg_bp_prob_cutoff,
                          hard_threshold=hard_threshold,
                          max_num_edges=max_num_edges,
                          no_lonely_bps=True,
                          nesting=True)
    fold_vectorize = make_fold_vectorize(complexity=3, nbits=15, fold=fold)
    from RaSE import compute_SNP_stability
    score = compute_SNP_stability(wild_seq, snp_tag=tag_tup, fold_vectorize=fold_vectorize)
    # print 'Rase-Score:', score
    return score


if __name__ == "__main__":
   main(sys.argv[1:])
