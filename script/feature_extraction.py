from __future__ import absolute_import
import argparse
import os
import sys
import re
import h5py
import glob
import multiprocessing
from scipy import interpolate
from statsmodels import robust
import traceback
import numpy as np
from tqdm import tqdm


def get_base_quality(reference,sam):
    """
    This function parse base qualities from sam file.

    Args:
        reference (str): The path to reference transcripts, in fasta format.
        sam (str): Alignment results from minimap2.

    Returns:
        dict: A dictionary containing base quality and sequence info.
    """
    base_quality_dict=dict()
    
    reference_dict=dict()
    with open(reference) as f:
        for line in f:
            line=line.rstrip()
            if ">" in line:
                
                contig=line.split()[0][1:]
                reference_dict[contig]=""
            else:
                sequence=line
                reference_dict[contig]=reference_dict[contig]+sequence

    with open(sam) as f:
        for line in f:
    
            if line[0] != "@":
                items =line.split("\t")
                id=items[0]
                flag=items[1]
                chr =items[2]
                start=int(items[3])
                CIGAR=items[5]
                seq=items[9]
                base_quality_string=items[10]
                base_quality_list=[ord(char)-33 for char in base_quality_string]   
                #Convert ascii string to int list. 
                mapped_base_quality_list=[]
                
                if  chr != "*":
                    
                    temp=""
                    index=0
                    new_seq=""
                    for char in CIGAR:
                        temp+=char
                        if char=="H":
                            num=int(temp[:-1])
                            temp=""
                        elif char=="S":
                            num=int(temp[:-1])
                            index+=num
                            temp=""
                        elif char=="M":
                            num=int(temp[:-1])
                            new_seq+=seq[index:index+num]
                            mapped_base_quality_list.extend(base_quality_list[index:index+num])
                            index+=num
                            temp=""
                        elif char=="I":
                            num=int(temp[:-1])
                            index+=num
                            temp=""
                        elif char=="D":
                            num=int(temp[:-1])
                            new_seq+="N"*num
                            mapped_base_quality_list.extend([0]*num)  
                            #padding 0 for deletions, as placeholder
                            temp=""

                    
                    if flag=="0":
                    #DRS reads only map to sense strand.
                        base_quality_dict[id] = [chr,start,reference_dict[chr][start-1:start+len(new_seq)-1],"|".join([str(x) for x in mapped_base_quality_list])]
                    

    return base_quality_dict


def get_events(fast5_path, basecall_group, basecall_subgroup,reverse = False):
	"""
    This function extract events from fast5 file.

    Args:
        fast5_path (str): The path to fast5 file.
        basecall_group (str): The default group from tombo output is "RawGenomeCorrected_000".
    	basecall_subgroup (str): The default basecall subgroup is "BaseCalled_template".
    Returns:
        dict: A dictionary containing base quality and sequence info.  to be done.
    """

	try:
		fast5_data = h5py.File(fast5_path, 'r')
	except IOError:
		raise IOError('Error opening file. Likely a corrupted file.')

	# Get raw data
	try:
		raw_data = list(fast5_data['/Raw/Reads/'].values())[0]
		
		raw_data = raw_data['Signal'][()]

		# ~ .value
	except Exception as e:
		
		raise RuntimeError(
			'Raw data is not stored in Raw/Reads/Read_[read#] so ' +
			'new segments cannot be identified.')

	# Read corrected data
	try:

		corr_data = fast5_data['/Analyses/'+basecall_group +'/' + basecall_subgroup + '/Events']
		corr_attrs = dict(list(corr_data.attrs.items()))
		corr_data = corr_data[()]

	except Exception as e:
        
		raise RuntimeError(('Corrected data not found.'))

	# Reading extra information
	corr_start_rel_to_raw = corr_attrs['read_start_rel_to_raw']  #
	if len(raw_data) > 99999999:
		raise ValueError(fast5_fn + ": max signal length exceed 99999999")
	if any(len(vals) <= 1 for vals in (corr_data, raw_data)):
		raise NotImplementedError(('One or no segments or signal present in read.'))
	event_starts = corr_data['start'] + corr_start_rel_to_raw
	event_lengths = corr_data['length']
	event_bases = corr_data['base']
	fast5_data.close()

	return raw_data, event_bases, event_starts, event_lengths

    
def get_signal(fast5_path):
	"""
	This function extracts the signal and sequence from a FAST5 file.

	Parameters:
	- fast5_path: The path to the FAST5 file.

	Returns:
	- line: A string containing the extracted information in a specific format.
	"""
	try:
		signal, sequence, signal_start, signal_length  = get_events(fast5_path, args.basecall_group,args.basecall_subgroup)
	except Exception as e:
		#print(152,e)
		return False, (None, None)
        
	signal = signal[::-1]
	signal = signal.tolist()
	signal_list=[]
	sequence="".join([x.decode() for x in sequence]) 

	for i in range(len(signal_length)):
		signal_list.append("*".join([str(x) for x in signal[signal_start[i]:signal_start[i]+signal_length[i]]]))
	
	line="%s\t%s\t%s\n"%(str(fast5_path).split("/")[-1].split(".")[0],sequence,"|".join(signal_list))     
    
	return line

def extract_feature(file_list):
	"""
	This function extracts signals from a list of files.

	Parameters:
	- file_list: A list of file paths.

	Returns:
	None
	"""
	base_quality_dict=get_base_quality(args.reference,args.sam)

	
	if int(args.process) > 1:
		results=[]
		pool = multiprocessing.Pool(processes = int(args.process))

		for file in file_list:

			result=pool.apply_async(get_signal,(file,))
			results.append(result)
		pool.close()

		pbar=tqdm(total=len(file_list),position=0, leave=True)
		nums=[]
		for result in results:
			num=result.get()

			if num:
				read_id,sequence,signal_string=num.split("\t")
				chr=base_quality_dict[read_id][0]
				start=int(base_quality_dict[read_id][1])
				reference_sequence=base_quality_dict[read_id][2]                        
				base_quality_list=base_quality_dict[read_id][3].split("|")

				feature=extract_5mer_features(read_id, chr, start, reference_sequence, base_quality_list, sequence, signal_string)
				process_id = os.getpid()
			pbar.update(1)

		pool.join()

	else:
		out=open(args.output,"w")
		nums=[]
		for file in file_list:

			num=get_signal(file)
			if num[0]:
				try:

					read_id,sequence,signal_string=num.split("\t")
					chr=base_quality_dict[read_id][0]
					start=int(base_quality_dict[read_id][1])
					reference_sequence=base_quality_dict[read_id][2]                        
					base_quality_list=base_quality_dict[read_id][3].split("|")
				
					feature=extract_5mer_features(read_id, chr, start, reference_sequence, base_quality_list, sequence, signal_string)
					
					if feature:
						out.write(feature)
				except Exception as e:
					print(traceback.print_exc())
					pass
		out.close()

 
	


def get_file_list(fast5_dir):
    cmd="find %s -name '*.fast5' >%s.txt" %(fast5_dir,fast5_dir)
    os.system(cmd)



def test1(args):
	print("extract")
	print(args.p)


def interp(x):
    """
    Interpolate and resample the current signal to make each base correspond to the same length of current signal.

    Args:
        x (float list): Origin current signal.

    Returns:
        float list: Reverse complement sequence.
    """
    l=len(x)
    y=x
    x=np.linspace(0,l-1,l)
    f=interpolate.interp1d(x,y,kind='slinear')
    x_new=np.linspace(0,l-1,100)
    y_new=f(x_new)
    y_new=np.around(y_new,4)
    return y_new.tolist()


def convert_base_name(base_name):
    """
    Converts a base name into a regular expression pattern.

    Args:
        base_name (str): Input base name to be converted.

    Returns:
        str: Regular expression pattern representing the converted base name.
    """
    merge_bases = {
        'A': 'A',
        'C': 'C',
        'G': 'G',
        'T': 'T',
        'M': '[AC]',
        'V': '[ACG]',
        'R': '[AG]',
        'H': '[ACT]',
        'W': '[AT]',
        'D': '[AGT]',
        'S': '[CG]',
        'B': '[CGT]',
        'Y': '[CT]',
        'N': '[ACGT]',
        'K': '[GT]'
    }
    pattern = ''
    for base in base_name:
        pattern += merge_bases.get(base, base)
    return pattern


def extract_5mer_features(read_id, chr, start, reference_sequence, base_quality_list, sequence, signal_string):
	"""
	Extracts features from a signal file.
    
	Args:
        signal_file (str): Path to the signal file.
        args (Namespace): Command-line arguments.
	"""

	kmer_fillter = convert_base_name(args.motif)
	site_list=""
	clip=10 
	count=0
	scaling="median_mad"
	try:


		raw_signal=[base_signal_string.split("*") for base_signal_string in signal_string.split("|")]   #string list
                            
		full_length_signal=np.array([x for x in re.split('\*|\|',signal_string)],dtype=int)   #for normlization
		full_length_signal_min=min(full_length_signal)
		full_length_signal_max=max(full_length_signal)
		full_length_signal_uniq=np.unique(full_length_signal)
                
		full_length_mean=np.mean(full_length_signal)
		full_length_std=np.std(full_length_signal)
                
		for index in range(clip,len(sequence)-clip):
                    
			kmer_sequence=sequence[index-2:index+3]        
			if len([x.start() for x in re.finditer(kmer_fillter,kmer_sequence)])==0:
				
				continue
			if sequence[index-2:index+3] != reference_sequence[index-2:index+3]:
                #print("sequence inconsitant error",read_id,sequence[index-2:index+3],reference_sequence[index-2:index+3])
				continue

			kmer_raw_signal=raw_signal[index-2:index+3]
			kmer_raw_signal=[np.array(x,dtype=float) for x in kmer_raw_signal]
                    
			if scaling=="min_max":
				kmer_raw_signal=[(x-full_length_signal_min)/(full_length_signal_max-full_length_signal_min) for x in kmer_raw_signal]     #min_max scaling
			elif scaling=="zscore":
				kmer_raw_signal=[(x-full_length_mean)/full_length_std for x in kmer_raw_signal]     #z-score scaling
			elif scaling=="median_mad":
				kmer_raw_signal=[(x-np.median(full_length_signal_uniq))/robust.mad(full_length_signal_uniq) for x in kmer_raw_signal]     #median_mad scaling
                    
			mean=[np.round(np.mean(x),3) for x in kmer_raw_signal]
			std=[np.round(np.std(x),3) for x in kmer_raw_signal]
			median=[np.round(np.median(x),3) for x in kmer_raw_signal]
			length=[len(x) for x in kmer_raw_signal]
			kmer_base_quality=base_quality_list[index-2:index+3]
                    
			for i in range(5):
				kmer_raw_signal[i]=interp(kmer_raw_signal[i])
                        
			line="%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(read_id,
                                                                                     chr,
                                                                                     start+index,
                                                                                     kmer_sequence,
                                                                                     "|".join([str(x) for x in  mean]),
                                                                                     "|".join([str(x) for x in std]),
                                                                                     "|".join([str(x) for x in  median]),
                                                                                     "|".join([str(x) for x in length]),
                                                                                     "|".join([str(x) for x in  kmer_base_quality]),
                                                                                     "|".join([str(x) for x in kmer_raw_signal[0]]),
                                                                                     "|".join([str(x) for x in kmer_raw_signal[1]]),
                                                                                     "|".join([str(x) for x in kmer_raw_signal[2]]),
                                                                                     "|".join([str(x) for x in kmer_raw_signal[3]]),
                                                                                     "|".join([str(x) for x in kmer_raw_signal[4]])
                                                                                    )
			site_list+=line
			count+=1
	except Exception as e:
		#print(traceback.print_exc())
		return 0
	return site_list

def main():
    
    fast5_files= glob.glob(os.path.join(args.input, '**/*.fast5'), recursive=True)

    extract_feature(fast5_files)

if __name__ == "__main__":
	#Create the main parse
	parser = argparse.ArgumentParser(description='Extract current signal from fast5 files.')
    


	parser.add_argument('-i', '--input', required = True, help="The directory containing fast5 files.")
	parser.add_argument('--basecall_group',default = "RawGenomeCorrected_000", help='The attribute group to extract the training data from. e.g. RawGenomeCorrected_000.')
	parser.add_argument('--basecall_subgroup', default='BaseCalled_template', help='Basecall subgroup Nanoraw resquiggle into. Default is BaseCalled_template.')
	parser.add_argument('-p','--process', default=1,help='Process.')
	parser.add_argument('--clip', default=10,help='The number of bases to be discarded at both ends.')
	parser.add_argument('-r','--reference',required = True,help='Reference transcripts fasta file.')
	parser.add_argument('--sam',required = True,help='Sam file.')
	parser.add_argument('-o', '--output', required = True, help="Output file.")
	parser.add_argument('-m', '--motif', required = True, help="Motif pattern to extract.")



	args = parser.parse_args()
    


	main()








